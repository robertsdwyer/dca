from msa import Msa
from dca import Dca
import numpy as np
import unittest

class TestDcaMethods(unittest.TestCase):
    
    def test_get_discount(self):
        msa = Msa('5296874.full.fas')
        msa.msa = msa.msa[0:3, :]
        self.assertEqual(len(Dca.get_discount(msa, 3, 0.9)), 3)
        self.assertEqual(Dca.get_discount(msa, 3, 0.9)[0], 1)
        self.assertEqual(Dca.get_discount(msa, 3, 0.9)[1], 0.5)
        
    def test_get_fiA_list(self):
        msa = Msa('5296874.full.fas')
        msa.msa =  msa.msa[0:3, :]
        discount_list = [1, 1, 1]
        alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N',
            'P','Q','R','S','T','V','W','Y','-']
        m = 3
        l = np.shape(msa.msa)[1]
        q = 21
        meff = sum(discount_list)
        fiA_list = Dca.get_fiA_list(msa, discount_list, alphabet, m, l, q, meff, 0)
        #test function w/o pseudocounts or discounts
        self.assertEqual(len(fiA_list), l*q)
        self.assertEqual(fiA_list[0], 0)
        self.assertEqual(fiA_list[5], 1/3.)
        self.assertEqual(fiA_list[20], 2/3.)
        self.assertEqual(fiA_list[(q*(l - 1))], 0)
        self.assertEqual(fiA_list[(q*(l - 1)) + 14], 1/3.)
        self.assertEqual(fiA_list[(q*(l - 1)) + 8], 2/3.)
        #test function w/ pseducounts but no discounts
        fiA_list = Dca.get_fiA_list(msa, discount_list, alphabet, m, l, q, meff, 0.5)
        self.assertEqual(len(fiA_list), l*q)
        self.assertEqual(fiA_list[0], 0.5*1./q)
        self.assertEqual(fiA_list[5], 0.5*1./q + 0.5*1/3.)
        self.assertEqual(fiA_list[20], 0.5*1./q + 0.5*2/3.)
        self.assertEqual(fiA_list[(q*(l - 1))], 0.5*1./q)
        self.assertEqual(fiA_list[(q*(l - 1)) + 14], 0.5*1./q + 0.5*1/3.)
        self.assertEqual(fiA_list[(q*(l - 1)) + 8], 0.5*1./q + 0.5*2/3.)
        #test function w/ discounts but no pseudocounts
        discount_list = [0.5, 1, 1]
        meff = 2.5
        fiA_list = Dca.get_fiA_list(msa, discount_list, alphabet, m, l, q, meff, 0)
        self.assertEqual(len(fiA_list), l*q)
        self.assertEqual(fiA_list[0], 0)
        self.assertEqual(fiA_list[5], 0.2)
        self.assertEqual(fiA_list[20], 0.8)
        self.assertEqual(fiA_list[(q*(l - 1))], 0)
        self.assertEqual(fiA_list[(q*(l - 1)) + 14], 0.2)
        self.assertEqual(fiA_list[(q*(l - 1)) + 8], 0.8)
        
    def test_get_fijAB_list(self):
        msa = Msa('5296874.full.fas')
        msa.msa =  msa.msa[0:3, 3:5]
        discount_list = [1, 1, 1]
        alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N',
            'P','Q','R','S','T','V','W','Y','-']
        m = 3
        l = np.shape(msa.msa)[1]
        q = 21
        meff = sum(discount_list)
        #test function w/o pseudocounts or discounts
        fijAB = Dca.get_fijAB_list(msa, discount_list, alphabet, l, m, q, meff, 0)
        self.assertEqual(np.shape(fijAB)[0], l*q - l)
        self.assertEqual(np.shape(fijAB)[1], l*q - l)
        self.assertEqual(fijAB[0, 0], 0)
        self.assertEqual(fijAB[2, 30], 1/3.)
        self.assertEqual(fijAB[30, 2], 1/3.)
        self.assertEqual(fijAB[14, 30], 1/3.)
        self.assertEqual(fijAB[30, 14], 1/3.)
        self.assertEqual(fijAB[3, 30], 1/3.)
        self.assertEqual(fijAB[30, 3], 1/3.)
        #test function w/ pseudocounts but no discounts
        fijAB = Dca.get_fijAB_list(msa, discount_list, alphabet, l, m, q, meff, 0.5)
        self.assertEqual(np.shape(fijAB)[0], l*q - l)
        self.assertEqual(np.shape(fijAB)[1], l*q - l)
        self.assertEqual(fijAB[2, 30], 0.5*1./(q*q) + 0.5*1/3.)
        self.assertEqual(fijAB[30, 2], 0.5*1./(q*q) + 0.5*1/3.)
        self.assertEqual(fijAB[14, 30], 0.5*1./(q*q) + 0.5*1/3.)
        self.assertEqual(fijAB[30, 14], 0.5*1./(q*q) + 0.5*1/3.)
        self.assertEqual(fijAB[3, 30], 0.5*1./(q*q) + 0.5*1/3.)
        self.assertEqual(fijAB[30, 3], 0.5*1./(q*q) + 0.5*1/3.)
        #test function w/o pseudocounts but with discounts
        discount_list = [0.5, 1.0, 1.0]
        meff = 2.5
        fijAB = Dca.get_fijAB_list(msa, discount_list, alphabet, l, m, q, meff, 0)
        self.assertEqual(np.shape(fijAB)[0], l*q - l)
        self.assertEqual(np.shape(fijAB)[1], l*q - l)
        self.assertEqual(fijAB[0, 0], 0)
        self.assertEqual(fijAB[2, 30], 0.2)
        self.assertEqual(fijAB[30, 2], 0.2)
        self.assertEqual(fijAB[14, 30], 0.4)
        self.assertEqual(fijAB[30, 14], 0.4)
        self.assertEqual(fijAB[3, 30], 0.4)
        self.assertEqual(fijAB[30, 3], 0.4)
        
    def test_get_cov_mat(self):
        fijAB = np.matrix([[1, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 0], [1, 0, 1, 0]])
        fiA = np.array([0.5, 0.5, 1, 0.5, 0.5, 1])
        alphabet = ['a', 'b', 'c']
        l = 2
        q = 3
        cov_mat = Dca.get_cov_mat(fiA, fijAB, alphabet, l, q)
        self.assertEqual(cov_mat[0, 0], 0.25)
        self.assertEqual(cov_mat[1, 1], 0.25)
        self.assertEqual(cov_mat[2, 2], 0.25)
        self.assertEqual(cov_mat[3, 3], 0.25)
        self.assertEqual(cov_mat[0, 1], 0)
        self.assertEqual(cov_mat[1, 0], 0)
        self.assertEqual(cov_mat[2, 3], 0)
        self.assertEqual(cov_mat[3, 2], 0)
        self.assertEqual(cov_mat[0, 2], 0.75)
        self.assertEqual(cov_mat[0, 3], 0.75)
        self.assertEqual(cov_mat[1, 2], 0.75)
        self.assertEqual(cov_mat[1, 3], 0.75)
        self.assertEqual(cov_mat[2, 0], 0.75)
        self.assertEqual(cov_mat[3, 0], 0.75)
        self.assertEqual(cov_mat[2, 1], 0.75)
        self.assertEqual(cov_mat[3, 1], -0.25)
        
    def test_get_submatrix(self):
        mtrx = np.matrix([[1, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 0], [1, 0, 1, 0]])
        alphabet = ['a', 'b', 'c']
        i = 0
        j = 0
        q = 3
        #check that first block is ok
        submtrx = Dca.get_submatrix(mtrx, alphabet, i, j, q)
        self.assertEqual(np.shape(submtrx)[0], 2)
        self.assertEqual(np.shape(submtrx)[1], 2)
        self.assertEqual(submtrx[0, 0], 1)
        self.assertEqual(submtrx[0, 1], 1)
        self.assertEqual(submtrx[1, 0], 1)
        self.assertEqual(submtrx[1, 1], 0)
        #check that last block is ok
        i = 1
        j = 1
        submtrx = Dca.get_submatrix(mtrx, alphabet, i, j, q)
        self.assertEqual(np.shape(submtrx)[0], 2)
        self.assertEqual(np.shape(submtrx)[1], 2)
        self.assertEqual(submtrx[0, 0], 0)
        self.assertEqual(submtrx[0, 1], 0)
        self.assertEqual(submtrx[1, 0], 1)
        self.assertEqual(submtrx[1, 1], 0)
        #check that mid block is ok
        i = 1
        j = 0
        submtrx = Dca.get_submatrix(mtrx, alphabet, i, j, q)
        self.assertEqual(np.shape(submtrx)[0], 2)
        self.assertEqual(np.shape(submtrx)[1], 2)
        self.assertEqual(submtrx[0, 0], 1)
        self.assertEqual(submtrx[0, 1], 1)
        self.assertEqual(submtrx[1, 0], 1)
        self.assertEqual(submtrx[1, 1], 0)
        
    def test_get_mean_field_weights(self):
        submtrx = np.matrix([[0, 1], [1, 0]])
        q = 3
        mfw = Dca.get_mean_field_weights(submtrx, q)
        self.assertEqual(np.shape(mfw)[0], 3)
        self.assertEqual(np.shape(mfw)[1], 3)
        self.assertEqual(mfw[0, 0], 1)
        self.assertEqual(mfw[0, 1], np.exp(-1))
        self.assertEqual(mfw[0, 2], 1)
        self.assertEqual(mfw[1, 0], np.exp(-1))
        self.assertEqual(mfw[1, 1], 1)
        self.assertEqual(mfw[1, 2], 1)
        self.assertEqual(mfw[2, 0], 1)
        self.assertEqual(mfw[2, 1], 1)
        self.assertEqual(mfw[2, 2], 1)
        
    def test_get_di_block(self):
        mfw = np.matrix([[0, 1], [1, 0]])
        fiA = np.array([1, 0.75, 0.5, 0.25])
        alphabet = ['a', 'b']
        mu1 = np.array([2, 3])
        mu2 = np.array([4, 5])
        i = 0
        j = 1
        q = 2
        fudge = 10**(-100)
        pdir = np.matrix([[12/45., 15/45.], [8/45., 10/45.]])
        pfac = np.matrix([[0.5, 0.25], [0.375, 0.75/4]])
        di_result = 0
        for k in range(2):
            for l in range(2):
                di_result += pdir[k, l]*np.log((pdir[k, l] + fudge)/(pfac[k, l] + fudge))
        self.assertEqual(Dca.get_di_block(mfw, fiA, alphabet, mu1, mu2, i, j, q), di_result)
        
    def test_map_key(self):
        alphabet = ['A', 'B', 'C']
        self.assertEqual(Dca.map_key(alphabet, 0, 'A', ref = True), 0)
        self.assertEqual(Dca.map_key(alphabet, 0, 'B', ref = True), 1)
        self.assertEqual(Dca.map_key(alphabet, 0, 'C', ref = True), 'skip')
        self.assertEqual(Dca.map_key(alphabet, 0, 'C', ref = False), 2)
        self.assertEqual(Dca.map_key(alphabet, 1, 'A', ref = True), 2)
        self.assertEqual(Dca.map_key(alphabet, 1, 'A', ref = False),3) 
        self.assertEqual(Dca.map_key(alphabet, 1, 'C', ref = True), 'skip')
        self.assertEqual(Dca.map_key(alphabet, 1, 'C', ref = False), 5)
        self.assertEqual(Dca.map_key(alphabet, 2, 'B', ref = True), 5)
        self.assertEqual(Dca.map_key(alphabet, 2, 'B', ref = False), 7)
    
    def test_get_fiA_drop_list(self):
        fiA = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
        alphabet = ['A', 'B', 'C']
        q = 3
        l = 3
        self.assertEqual(len(Dca.get_fiA_drop_list(fiA, alphabet, l, q)), 6)
        self.assertEqual(Dca.get_fiA_drop_list(fiA, alphabet, l, q)[5], 8)
        
if __name__ == '__main__':
    unittest.main()