import os

class Sequence:
    
    class Fasta:
        def __init__(self,inp):
            
            self.inp = inp
            self.raw_fasta = self.parse_input()
            self.parse_fasta()
            
        def write_fasta(self, filename):
            #check whether filename exists
            exists = os.path.exists(filename)
            if exists:
                raise ValueError(filename + ' already exists')
            else:
                #create fasta string
                fasta_string = self.seq_to_string()
                #write fasta
                f = open(filename, 'w')
                f.write(fasta_string)
                f.close()

        def seq_to_string(self):
            #determine how many iterations to loop for linebreaks
            iterations = self.length()/70
            if self.length()%70 > 0:
                iterations += 1
            #loop through sequence adding linebreaks
            fasta_chunked = []
            for i in range(iterations):
                seq_chunk = self.fasta.seq[(i*70):((i+1)*70)]
                fasta_chunked.append(seq_chunk)
            fasta_string = '\n'.join([self.fasta.header] + fasta_chunked)
            return fasta_string
                
        def parse_fasta(self):
            #define fasta alphabet
            fasta_alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N',
            'P','Q','R','S','T','V','W','Y','-']
            #check whether raw_fasta is an object of type str
            if type(self.raw_fasta) == str:
                #check whether raw_fasta has header marker >
                if '>' in self.raw_fasta:
                    #split raw_fasta by lines
                    fasta_lines = self.raw_fasta.split('\n')
                    #define header as first line
                    self.header = fasta_lines[0]
                    #check that first line (header) contains header marker >
                    if '>' not in self.header:
                        raise ValueError('Header mis-parsed: does not contain >'
                        )
                    #create sequence from remaining lines
                    self.seq = ''.join(fasta_lines[1:])
                    #convert all to uppercase
                    self.seq = self.seq.upper()
                    if ('X' in self.seq):
                        self.seq = self.seq.replace('X','-')
                    for aa in self.seq:
                        if aa not in fasta_alphabet:
                            raise ValueError(' '.join([
                            'Inadmissable amino acid code', aa, 'found in seq'
                            ]))
                            
        def parse_input(self):
            #where length of self.ipn is small, assume a filename reference
            if len(self.inp) < 100:
                fasta = self.get_nucleator_fasta_text()
                if not fasta:
                    fasta = self.inp
            else:
                #assume longer sequences are raw fasta files
                fasta = self.inp  
            return fasta
     
        def get_nucleator_fasta_text(self):
            #Checks whether self.inp refers to full file name, extension-less file-
            #name, or filename prefix, i.e. protein name, in lib/nucleators; if so 
            #loads and returns object of class Fasta; else returns False
            #
            #instantiate nucleator_filename
            nucleator_filename = False
            #set relative path to nucleator lib
            nucleator_path = 'lib/nucleators/'
            #obtain list of files in nucleator_path directory
            nucleator_files = os.listdir(nucleator_path)
            #check whether inp in nucleator_files
            if self.inp in nucleator_files:
                nucleator_filename = nucleator_files[nucleator_files.index(self.inp)
                ]
            else:
                #iterate through files and remove file extension
                nucleator_files_trunc = []
                for f in nucleator_files:
                    nucleator_files_trunc.append(f.split('.')[0])
                #check whether self.inp matches any truncated nucleator files
                if self.inp in nucleator_files_trunc:
                    #pull original path
                    nucleator_filename = nucleator_files[
                    nucleator_files_trunc.index(self.inp)]
                else:
                    #iterate through files and pull prefix, i.e. protein name
                    nucleator_files_prefix = []
                    for f in nucleator_files:
                        nucleator_files_prefix.append(f.split('_')[0])
                    #check whether self.inp in nucleator_files_prefix
                    if self.inp in nucleator_files_prefix:
                        #check whether prefix is unique
                        if nucleator_files_prefix.count(self.inp) > 1:
                            raise ValueError(
                            'Multiple files found with prefix {0}'.format(self.inp))
                        else:
                            nucleator_filename = nucleator_files[
                            nucleator_files_prefix.index(self.inp)]
            #check whether nucleator_filename has been defined
            if nucleator_filename:
                #load filename
                with open('/'.join([nucleator_path,nucleator_filename]),'r') as f:
                    fasta = f.read()
            else:
                fasta = False
            return fasta

    def __init__(self,inp):
        self.inp = inp
        self.fasta = self.Fasta(self.inp)
        
    def length(self):
        # returns the length of a sequence
        length = len(self.fasta.seq)
        return length
        
    def get_resi(self, idx, oneidx=False):
        #returns the identity of a residue at a given idx
        #
        #determine whether zero-indexed (default) or one-indexed
        if oneidx:
            idx = idx - 1
        return self.fasta.seq[idx]

            
            
            
            
                
            
            
        