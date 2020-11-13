##
# Testes para py
#
# Todos os testes no folder "test" são executados quando a ação "Test" é invocada.
#
##


import unittest
from sequencinator.intake import ler_FASTA_seq, ler_seq


class MyTestCase(unittest.TestCase):
    
    def test_ler_seq(self):

        f=open(r"seq_testes/seq_testar1.txt")
        self.assertEqual(ler_seq(f), "SVRNFLQRVVLGVDTKKKNSJDGIABJAHIUHUSNJNIYDDEINSLDTEQINDGAPVEKKSPLGYNVSDFTACYLIIEGVIGTGIFATPATILKSVGSVGASYVFWCVGSHUSHUJSUYSVBTSKSJSOJDOJSYFRRRSGAQVAYLEQAYPKPKFLVPVLYAAVSVTLSDDNSNDS")
        f.close()
        
        f=open(r"seq_testes/seq_testar2")
        self.assertEqual(ler_seq(f), "ATGCGTGAGCTAGGAAACGATAGACGAGTAGACCGATATGC")
        f.close()

        f=open(r"seq_testes/seq_testar3")
        self.assertEqual(ler_seq(f), "ATCGATCGATGACGATGCAGGTGCTAGTGCAGACTGTCTGCT")
        f.close()

        f=open(r"seq_testes/seq_testar4")
        self.assertEqual(ler_seq(f), "")
        f.close()

    def test_ler_FASTA_seq(self):
        
        f=open(r"FASTA_testes/FASTA_testar1.txt")
        self.assertEqual(ler_FASTA_seq(f), "MSSAFVSVRNFLQRVVLGVDTKKNIYDDEINSLDTEQINDGAPVEKKSPLGYNVSDFTACYLIIEGVIGTGIFATPATILKSVGSVGASYVFWCVGFVVNQFTVLMYVEYVTYFRRRSGAQVAYLEQAYPKPKFLVPVLYAAVSVTLSYITSSAASFSQYVFEGANYEATAWQQRGLSILPLFLAAIFTTLSTKWTLRLNSIIGWCKVCFIFFIAFSGFAALAGSTKAPKNHDIFKNAWEGTTTDGNSISNAILKVVFSFGGAPYAFTVVAETHPKNTIK")
        f.close()
       
        f=open(r"seq_testes/seq_testar2")
        self.assertEqual(ler_FASTA_seq(f), "ATGCGTGAGCTAGGAAACGATAGACGAGTAGACCGATATGC")
        f.close()
        
        f=open(r"FASTA_testes/FASTA_testar3")
        self.assertEqual(ler_FASTA_seq(f), "CGUCAGGUAGCUCAGUCGCGGACUGACUUACGCGUACGUCGAUAGCUGCUCGUGAC")
        f.close()

        f=open(r"FASTA_testes/FASTA_testar4")
        self.assertEqual(ler_FASTA_seq(f), "")
        f.close()


if __name__ == '__main__':
    unittest.main()

