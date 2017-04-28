from scipy.spatial.distance import hamming
import numpy as np

class Dca:
    def __init__(self, msa, similarity_discount = 0.75, lam = 0.5):
        #check that msa is of type Msa
        if msa.__class__.__name__ != 'Msa':
            raise ValueError('Argument msa must be of class Msa')
        self.msa = msa
        self.alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N',
            'P','Q','R','S','T','V','W','Y','-']
        self.m = self.msa.depth()
        self.l = self.msa.length()
        self.q = len(self.alphabet)
        self.similarity_discount = similarity_discount
        self.discount_list = self.get_discount(self.msa, self.m, self.similarity_discount)
        self.meff = sum(self.discount_list)
        self.lam = lam
        self.fiA_list = self.get_fiA_list(self.msa, self.discount_list, self.alphabet, self.m, self.l, self.q, self.meff, self.lam)
        fijAB = self.get_fijAB_list(self.msa, self.discount_list, self.alphabet,  self.m, self.l, self.q, self.meff, self.lam)
        cov_matrix = self.get_cov_mat(self.fiA, fijAB, self.alphabet, self.l, self.q)
        invc = np.linalg.inv(cov_matrix)
        self.di_list = self.get_di()
    @staticmethod
    def get_discount(msa, m, similarity_discount):
        discount_list = [1]*m
        if similarity_discount < 1:
            for i in range(m):
                seqi = msa.msa[i, :]
                for j in range(i + 1, m):
                    seqj = msa.msa[j, :]
                    ham = hamming(seqi, seqj)
                    if (1 - ham) > similarity_discount:
                        discount_list[i] += 1
                        discount_list[j] += 1
            discount_list = 1./np.array(discount_list)
        return discount_list
    @staticmethod   
    def get_fiA_list(msa, discount_list, alphabet, m, l, q, meff, lam):
        fiA_vect = np.array([0.]*q*l)
        for i in range(m):
            for j in range(l):
                indj = Dca.map_key(alphabet, j, msa.msa[i, j], ref = False)
                fiA_vect[indj] += discount_list[i]
        fiA_vect = fiA_vect/float(meff)
        fiA_vect = (1 - lam)*fiA_vect + lam/float(q)
        return fiA_vect
    @staticmethod    
    def get_fijAB_list(msa, discount_list, alphabet, l, m, q, meff, lam):
        dim = (q - 1)*l
        fijAB_matrix = np.zeros((dim, dim))
        for i in range(m):
            for j in range(l):
                indj = Dca.map_key(alphabet, j, msa.msa[i, j])
                if indj != 'skip':
                    for k in range(j + 1, l):
                        indk = Dca.map_key(alphabet, k, msa.msa[i, k])
                        if indk != 'skip':
                            fijAB_matrix[indj, indk] += float(discount_list[i])
        fijAB_matrix = fijAB_matrix/float(meff)
        fijAB_matrix = fijAB_matrix + np.transpose(fijAB_matrix)      
        fijAB_matrix = (1 - lam)*fijAB_matrix + lam/float(q*q)
        return fijAB_matrix
    @staticmethod
    def get_cov_mat(fiA, fijAB, alphabet, l, q):
        fiA_drop = Dca.get_fiA_drop_list(fiA, alphabet, l, q)
        c_matrix = fijAB - np.transpose(np.matrix(fiA_drop))*fiA_drop
        # convert diagonal blocks to zeros
        for i in range(l):
            start_block = Dca.map_key(alphabet, i, alphabet[0])
            end_block = start_block + q - 1
            c_matrix[start_block:end_block, start_block:end_block] = 0
        diag = []
        for i in range(l):
            start_slice = Dca.map_key(alphabet, i, alphabet[0], ref = False)
            end_slice = start_slice + q - 1
            diag += list(fiA[start_slice:end_slice])
        diag = np.array(diag)
        diag = diag*(1 - diag)
        np.fill_diagonal(c_matrix, diag)
        return c_matrix
    @staticmethod    
    def get_submatrix(invc, alphabet, i, j, q):
        start_blocki = Dca.map_key(alphabet, i, alphabet[0])
        end_blocki = start_blocki + q - 1
        start_blockj = Dca.map_key(alphabet, j, alphabet[0])
        end_blockj = start_blockj + q - 1
        submatrix = invc[start_blocki:end_blocki, start_blockj:end_blockj]
        return submatrix
    @staticmethod
    def get_mean_field_weights(submatrix, q):
        mfw = np.zeros((q, q)) + 1
        mfw[0:np.shape(submatrix)[0], 0:np.shape(submatrix)[1]] = np.exp(-submatrix)
        return mfw
    @staticmethod
    def compute_mu(mfw, alphabet, i, j, q, fiA):
        pi_start = Dca.map_key(alphabet, i, alphabet[0], ref = False)
        pi_end = pi_start + q
        pj_start = Dca.map_key(alphabet, j, alphabet[0], ref = False)
        pj_end = pj_start + q
        
        pi = fiA[pi_start, pi_end]
        pj = fiA[pj_start, pj_end]
        
        epsilon = 0.0001
        diff = 1.0
        mu1 = np.zeros((1, q))
        mu2 = np.zeros((1, q))
        while diff > epsilon:
            scra1 = mu2*np.transpose(mfw)
            scra2 = mu1*mfw
            new1 = pi/scra1
            new1 = pi/sum(pi)
            new2 = pj/scra2
            new2 = pj/sum(pj)
            diff = max(max(new1), max(new2))
            mu1 = new1
            mu2 = new2
        return mu1, mu2
    
    @staticmethod
    def get_di(invc, fiA, alphabet, i, j, l, q):
        di_list = []
        for i in range(l):
            for j in range(i+1, l):
                submatrix = Dca.get_submatrix(invc, alphabet, i, j, q)
                mfw = Dca.get_mean_field_weights(submatrix, q)
                mu1, mu2 = Dca.compute_mu(mfw, alphabet, i, j, q, fiA)
                di_ij = Dca.get_di_block(mfw, fiA, alphabet, mu1, mu2, i, j, q)
                di_list.append(tuple(di_ij, i, j))
        di_list = sorted(di_list, reverse = True)
        return di_list
    @staticmethod             
    def get_di_block(mfw, fiA, alphabet, mu1, mu2, i, j, q):
        pi_start = Dca.map_key(alphabet, i, alphabet[0], ref = False)
        pi_end = pi_start + q
        pj_start = Dca.map_key(alphabet, j, alphabet[0], ref = False)
        pj_end = pj_start + q
        pi = fiA[pi_start:pi_end]
        pj = fiA[pj_start:pj_end]
        
        mu1 = np.matrix(mu1)
        mu2 = np.matrix(mu2)
        
        fudge = 10**(-100)
        pdir = mfw * (np.transpose(mu1) * mu2)
        pdir = pdir/float(np.sum(pdir))
        pfac = np.transpose(np.matrix(pi)) * np.matrix(pj)
                
        di = sum( np.diagonal( np.transpose(pdir) * np.log(( pdir + fudge) / (pfac + fudge))))
        return di
    @staticmethod  
    def map_key(alphabet, i, val, ref = True):
        alphabet_len = len(alphabet)
        if ref:
            alphabet_len = alphabet_len - 1
        blockstart = i*alphabet_len
        ind = alphabet.index(val)
        if ind < alphabet_len:
            ind = blockstart + ind
        else:
            ind = 'skip'
        return ind
    @staticmethod
    def get_fiA_drop_list(fiA, alphabet, l, q):
        fiA_drop_list = []
        for i in range(l):
            pi_start = Dca.map_key(alphabet, i, alphabet[0], ref = False)
            pi_end = pi_start + q - 1
            pi = fiA[pi_start:pi_end]
            fiA_drop_list += list(pi)
        fiA_drop_list = np.array(fiA_drop_list)
        return fiA_drop_list