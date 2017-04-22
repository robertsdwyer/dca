import subprocess
import os
from sequence import Sequence

class Msa:
    
    def __init__(self,sequence):
        # check that sequence is of type Sequence
        if sequence.__class__.__name__ == 'Sequence':
            self.msa = self.get_msa()
        else:
            # check for known extension
            extension = sequence[-4:]
            if extension == '.fas':
                self.msa = self.load_msa(sequence)
                self.check_msa(self.msa)
            elif extension == '.a3m':
                self.msa = self.load_msa(sequence)
            else:
                raise ValueError("Msa class must be initiated with an existing MSA or a sequence of type 'Sequence'")
            
    def get_msa(self):
        # create new msa filename
        msa_path = '/Users/hannahhartig/Documents/dca/lib/msas'
        msa_filename = self.sequence.fasta.header.replace('<','').replace(' ','') + '.a3m'
        # check to see whether msa filename exists
        msa_filenames = os.listdir(msa_path)
        # instantiate getmsa as True
        getmsa = True
        if msa_filename in msa_filenames:
            #switch getmsa to False
            getmsa = False
            #get user input on whether to overwrite
            getting_input = True
            while getting_input:
                overwrite = raw_input('Overwrite existing msa ' + msa_filename + '? (y/n)')
                if overwrite == 'y':
                    getmsa = True
                    getting_input = False
                elif overwrite == 'n':
                    getmsa = False
                else:
                    print 'command ' + overwrite + ' not recognized; choose one of y/n'
        if getmsa:
            #create full msa_filename
            msa_filename = '/'.join([msa_path, msa_filename])
            #create temporary .fasta file
            fasta_temp_path = '/Users/hannahhartig/Documents/dca/lib/tmp'
            fasta_temp_filename = 'tmp.fasta'
            fasta_temp_filename_full = '/'.join([fasta_temp_path, fasta_temp_filename])
            #check to see whether tmp.fasta exists and if so, delete it
            fasta_temps = os.listdir(fasta_temp_path)
            if fasta_temp_filename in fasta_temps:
                os.remove(fasta_temp_filename_full)
            #save temporary .fasta file
            self.sequence.write_fasta(fasta_temp_filename_full)
            #write msa
            hhblits_command = 'hhblits -cpu 2 -i ' + fasta_temp_filename_full + ' -d databases/uniport20 -oa3m ' + msa_filename
            subprocess.check_output(hhblits_command)
            #load msa
            self.msa.a3m = self.load_msa(msa_filename)
            #convert to fas
            msa_filename_fas = msa_filename[0:(len(msa_filename)-3)] + 'fas'
            reformat_command = 'reformat.pl a3m fas ' + msa_filename + ' ' + msa_filename_fas
            subprocess.check_output(reformat_command)
            #load msa fas
            self.msa = self.load_msa(msa_filename_fas)
 
    @staticmethod   
    def load_msa(filename):
        # check whether filenam is str
        if not type(filename) == str:
            raise ValueError('Argument filename must be of type str')
        #check whether filename exists
        file_path = '/Users/hannahhartig/Documents/dca/lib/msas'
        if file_path not in filename:
            filename = '/'.join([file_path, filename])
        if not os.path.exists(filename):
            raise ValueError('File ' + filename + ' does not exist')
        else:
            with open(filename,'r') as f:
                    msa = f.read()
        return msa
        
    @staticmethod
    def check_msa(msa_string):
        #check that msa is an msa
        if not type(msa_string) == str:
            raise ValueError('Argument msa_string must be of type string')
        elif '>' not in msa_string:
            raise ValueError('Argument msa_string must be fasta format containing character >')
        elif msa_string.count('>') <= 1:
            raise ValueError('Argument msa_string must contain multiple sequences delimited by >')
        msa_split = msa_string.split('>')
        # check that all msa elements are valid sequences of the same length
        ln_last = 0
        for seq in msa_split:
            seq = '>' + seq
            seq = Sequence(seq)
            if ln_last == 0:
                ln_last = seq.length()
            else:
                if ln_last != seq.length():
                    raise ValueError('Argument msa_string must have sequences of equal lengths')   
        return True
        
    def dca(self):
        
        
        

