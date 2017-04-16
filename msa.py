import subprocess
import os

class Msa:
    
    def __init__(self,sequence):
        # check that sequence is of type Sequence
        if not sequence.__class__.__name__ == 'Sequence':
            raise ValueError("Msa class must be initiated with a sequence of type 'Sequence'")
        else:
            self.sequence = sequence
        self.msa = self.get_msa()
        
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
            self.msa.fas = self.load_msa(msa_filename_fas)
            
    def load_msa(self, filename):
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
