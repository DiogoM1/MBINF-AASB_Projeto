import unittest

from sequencinator.sequences import complemento_inverso, transcricao, traducao, reading_frames, valida, contar_bases
from sequencinator.protein import get_proteins


class teste_funcoes(unittest.TestCase):
    
    def test_ADN(self):        
        self.assertTrue(valida("ATGCTGCAGTGAT"))
        self.assertFalse(valida("AGTPAGTPP"))
    
    def test_RNA(self):
        self.assertFalse(valida("AUGGUUUCA"))
        
    def test_minusculas(self):
        self.assertTrue(valida("agttactacgattatg"))
        self.assertFalse(valida("agttacmknhtacaawgatdstaetg"))
    
    def test_mistura(self):
        self.assertTrue(valida("taAaCaaTccGg"))
        self.assertFalse(valida("AuGdUUAa_?/5Caa")) #testa para seq. invalida   
  
    def test_bases(self):
       self.assertEqual(contar_bases("ACGTAG"), {"A":2, "C":1, "G":2, "T":1}) #testa seq normal
       self.assertNotEqual(contar_bases("ACGTAG"), {"A":5, "C":3, "G":10, "T":45}) 
    
    def test_minusculas_b(self):   
       self.assertEqual(contar_bases("atgtttcagaaag"), {"A":5, "C":1, "G":3, "T":4})
        
    def test_nulo(self):        
        self.assertIsNone(contar_bases("AUGGUUUCA")) #testa para RNA
        self.assertIsNone(contar_bases("Affhth_reeet21?*")) #testa para seq. invalida
        
    def test_maimin(self):
        self.assertEqual(contar_bases("AgtaCGTcAG"), {"A":3, "C":2, "G":3, "T":2})
        self.assertEqual(contar_bases("agtaaCGTcgaG"), {"A":4, "C":2, "G":4, "T":2}) #maiusculas e minusculas
        
    def test_inverso(self):
        self.assertEqual(complemento_inverso("GTACTGATGCAAAACCGGTT"), "AACCGGTTTTGCATCAGTAC") #testa seq normal
        self.assertEqual(complemento_inverso("TCgcGCGCaTgggCC"), "GGCCCATGCGCGCGA") #maiusculas e minusculas
        self.assertEqual(complemento_inverso("tgctgtagtcagtcgtagtcata"), "TATGACTACGACTGACTACAGCA") #minusculas
        self.assertFalse(complemento_inverso("taUAaCauaTccGg")) #testa para RNA
        self.assertFalse(complemento_inverso("AuGdUUAa_?/5Caa")) #testa para seq. invalida
        
    def test_transcricao(self):
        self.assertEqual(transcricao("TCTATCTCCGGTCCAAGTCA"), "UCUAUCUCCGGUCCAAGUCA")
        self.assertEqual(transcricao("atcgatgatcgtagtctacggtcat"), "AUCGAUGAUCGUAGUCUACGGUCAU")
        self.assertEqual(transcricao("TcTaTCTCCggTCCAagTcA"), "UCUAUCUCCGGUCCAAGUCA")
        self.assertFalse(transcricao("aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa")) #testa para seq. invalida
        self.assertFalse(transcricao("taUAaCauaTccGg")) #testa para RNA
        
        
    def test_traducao(self):
        self.assertEqual(traducao("atgtagctgatcgtagctagtctatc"), "M_LIVASL")
        self.assertEqual(traducao("GCAATACCGTGGTCTAGCCTCAAATTCGGCTCTGTTACCTCGAGCGTTATGTGTCAAATGGCG"), "AIPWSSLKFGSVTSSVMCQMA")
        self.assertEqual(traducao("TcTaTCTCCAtgggTCCAagTcAttaA"), "SISMGPSH_") #maiusculas e minusculas
        self.assertFalse(traducao("aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa")) #testa para seq. invalida
        self.assertFalse(traducao("UauaaCuaGucgucuaaUACuccacucGg")) #testa para RNA
        
    def test_reading_frames(self):
        self.assertCountEqual(reading_frames("AGTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA"), ['AGTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'TGCGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT', 'GTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'GCGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT', 'TGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'CGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT'])
        self.assertCountEqual(reading_frames("tgcgggcagattatgtaggttgagagatgcgggagaagttctcgaccttcccgtgggacgtgaacctatcccctaatagagcattccgttcgagcatggc"), ['TGCGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'GCCATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA', 'GCGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'CCATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA', 'CGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'CATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA'])
        self.assertFalse(reading_frames("rueth8347y32894u59806y4957813r837gh25gh")) #testa para seq. invalida


    def test_get_proteins(self):
        self.assertEqual(get_proteins("ATGACCGTAA"), [])
        self.assertEqual(get_proteins("GCAATACCGTGGTCTAGCCTCAAATTCGGCTCTGTTACCTCGAGCGTTATGTGTCAAATGGCG"), ["AIPWSSLKFGSVTSSVMCQMA"])
        self.assertEqual(get_proteins("TcTaTCTCCAtgggTCCAagTcAttaA"), ["SISMGPSH_"])
        self.assertFalse(get_proteins("aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa")) #testa para seq. invalida
        self.assertFalse(get_proteins("UauaaCuaGucgucuaaUACuccacucGg")) #testa para RNA
        

       
        
if __name__ == "__main__":
    unittest.main()