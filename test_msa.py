import unittest
from msa import Msa

class TestMsaMethods(unittest.TestCase):
    
    def test_load_msa(self):
        msa_filename_fileinitfas = '5296874.full.fas'
        #check that msa is loaded properly from filename
        self.assertEqual(Msa.load_msa(msa_filename_fileinitfas)[0:5], '>4HTS')
        #check that load_msa throws exceptions for invalid inputs
        with self.assertRaises(ValueError):
            Msa.load_msa(3)
            Msa.load_msa('notafilename')
            Msa.load_msa('notafilename.fas')
            
    def test_check_msa(self):
        msa_filename_fileinitfas = '5296874.full.fas'
        msa_example = '''>4HTS:A|PDBID|CHAIN|SEQUENCE
GSVDMPL--TE-----HLREL----RYRLIISIIA------FLI--GS-GIA---------F----Y---------F-----------------A-K---
-----------Y-----------V-F-----------------------E---------I--------------------------L-------------
A------AP---------------------I----------------------L-----------K----------------SY----------------
-----------------------------P----------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
------------------------------------------------E---------------------------------------------------
---------V-----------------E---------L------------------IT------------------LSPT-------EP-----------
---LF-ILIK--ISLAVGFIIAS-PVILYQF----------WRFIE-PA-----LY-----SH--E----KRAFIP-----LL----------LGSI-LL
-FMLGALFAYFIVLPLALKFLL------------------------------------------------------------------------------
------------------------------------------------G-LG-------F---------T------------------------------
--------Q-------------------------------L----------------------------L---A---T----------------------
-------------P------------------------------------------------------------------------------------Y-
-----L-------SV---DM------YISF--------VLKL-VVAFG----IAFE-MPIVLY-V-LQKAGVIT-PEQL-A-------------------
-------S---FR-KYFIV-IAFVIG-A-I-I--A--P--DV--S--TQ------VLMAI--P-LLLLYE----ISIFLGK--LATR
>tr|C8X1I2|C8X1I2_DESRD Sec-independent protein translocase protein TatC OS=Desulfohalobium retbaense (strain DSM 5692) GN=tatC PE=3 SV=1
---RMGL--ID-----HLNDL----RKSIVRSVVA------ALV--GM-VAC---------Y----A---------F-----------------A-Q---
-----------K-----------L-F-----------------------D---------Y--------------------------L-------------
M------LP---------------------L----------------------Y-----------K----------------AL----------------
-----------------------------Pe---------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
------------------------------------------------G---------------------------------------------------
---------S-----------------T---------L------------------IY------------------TAPH-------EA-----------
---FF-TYIK--VAFVAGLFLTS-PFIFYQF----------WSFVA-PG-----LY-----KH--E----RKWLIP-----IA----------FFSA-AF
-FCTGAIFGYSIVFPWGYKFFM------------------------------------------------------------------------------
------------------------------------------------G----------F---------A------------------------------
--------E-------------------------------D----------------------------L---I---R----------------------
-------------P------------------------------------------------------------------------------------M-
-----L-------TM---RE------AFSF--------AMRL-LIAFG----VVFE-LPLVIF-F-LARLGLVN-AAWL-R-------------------
-------K---KR-KYAIL-IAFILS-A-L-L--T--Pp-DM--V--TQ------SFMAG--P-LALLYE----LSIWVAA--VFGK'''
        self.assertTrue(Msa.check_msa(msa_example))
        msa = Msa.load_msa(msa_filename_fileinitfas)
        self.assertTrue(Msa.check_msa(msa))
        with self.assertRaises(ValueError):
            Msa.check_msa(3)
            Msa.check_msa(msa_example.split('>')[0])
            Msa.check_msa('>' + '/n'.join(msa_example.split('>')))
            Msa.check_msa('>' + '/n>'.join([msa_example.split('>')[0], msa_example.split('>')[1][0:-1]]))
            
    def test_msa(self):
        msa_filename_fileinitfas = '5296874.full.fas'
        msa = Msa(msa_filename_fileinitfas)
        #check that msa is loaded properly from filename
        self.assertEqual(msa.msa[0:5], '>4HTS')
        
if __name__ == '__main__':
    unittest.main()