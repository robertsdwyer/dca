import unittest
from sequence import Sequence

class TestSeqMethods(unittest.TestCase):
    
    def test_fasta(self):
        fasta1 = '''>NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]
MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQAD
TNQEQAKTEDAKRHDKEQV'''
        #fasta should fail because there is an aa X at position 8
        fastafail1 = '''>NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]
MGGISIWXLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQAD
TNQEQAKTEDAKRHDKEQV'''
        #fasta should fail because there is no > to start the header
        fastafail2 = '''NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]
MGGISIWXLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQAD
TNQEQAKTEDAKRHDKEQV'''
        #fasta should fail because there is an additional > outside the header
        fastafail3 = '''>NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]
>MGGISIWXLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQAD
TNQEQAKTEDAKRHDKEQV'''
        #Check that Fasta parses data properly
        self.assertEqual(Sequence.Fasta(fasta1).header, '>NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]')
        self.assertEqual(Sequence.Fasta(fasta1).seq, 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV')
        self.assertEqual(Sequence.Fasta('TatA_Ecoli_K12_MG1655.txt').seq, 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV')
        self.assertEqual(Sequence.Fasta('TatA_Ecoli_K12_MG1655').seq, 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV')
        self.assertEqual(Sequence.Fasta('TatA').seq, 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV')
        #Check that Fasta raises appropriate errors
        with self.assertRaises(ValueError):
            Sequence.Fasta(fastafail1)
            Sequence.Fasta(fastafail2)
            Sequence.Fasta(fastafail3)
            
    def test_seq(self):
        fasta1 = '''>NP_418280.4 TatABCE protein translocation system subunit [Escherichia coli str. K-12 substr. MG1655]
MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQAD
TNQEQAKTEDAKRHDKEQV'''
        self.assertEqual(Sequence(fasta1).fasta.seq, 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV')
        self.assertEqual(Sequence(fasta1).length(), 89)
        self.assertEqual(Sequence(fasta1).get_resi(0), 'M')
        self.assertEqual(Sequence(fasta1).get_resi(88), 'V')
        self.assertEqual(Sequence(fasta1).get_resi(7), 'Q')
        self.assertEqual(Sequence(fasta1).get_resi(1, oneidx = True), 'M')
        self.assertEqual(Sequence(fasta1).get_resi(89, oneidx = True), 'V')
        self.assertEqual(Sequence(fasta1).get_resi(8 ,oneidx = True), 'Q')
          
if __name__ == '__main__':
    unittest.main()