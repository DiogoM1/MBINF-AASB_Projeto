import unittest

from sequencinator.sequences import complemento_inverso, transcricao, traducao, reading_frames, valida, contar_bases, valida_rna
from sequencinator.protein import get_proteins


class teste_funcoes(unittest.TestCase):
    
    def test_valida(self):
        self.assertTrue(valida("ATGCTGCAGTGAT")) #ADN
        self.assertFalse(valida("AGTPAGTPP")) #ADN
        self.assertFalse(valida("AUGGUUUCA")) #RNA
        self.assertTrue(valida("agttactacgattatg"))  #minusculas
        self.assertFalse(valida("agttacmknhtacaawgatdstaetg"))  #minusculas
        self.assertTrue(valida("taAaCaaTccGg")) #maiusculas e minusculas
        self.assertFalse(valida("AuGdUUAa_?/5Caa")) #testa para seq. invalida

    def test_valida_RNA(self):
        self.assertFalse(valida_rna("ATGCTGCAGTGAT")) #ADN
        self.assertFalse(valida_rna("AGTPAGTPP")) #ADN
        self.assertTrue(valida_rna("AUGGUUUCA")) #RNA
        self.assertTrue(valida_rna("aguuacuacgauuaug"))  #minusculas
        self.assertFalse(valida_rna("agttacmknhtacaawgatdstaetg"))  #minusculas
        self.assertTrue(valida_rna("UaAaCaaUccGg")) #maiusculas e minusculas
        self.assertFalse(valida_rna("AuGdUUAa_?/5Caa")) #testa para seq. invalida
  
    def test_contar_base(self):
        self.assertEqual(contar_bases("ACGTAG"), {"A":2, "C":1, "G":2, "T":1}) #testa seq normal
        self.assertNotEqual(contar_bases("ACGTAG"), {"A":5, "C":3, "G":10, "T":45})
        self.assertEqual(contar_bases("atgtttcagaaag"), {"A":5, "C":1, "G":3, "T":4})  #minusculas
        self.assertEqual(contar_bases ("AUGGUUUCA"), {"A":2, "C":1, "G":2, "U":4}) #testa para RNA
        self.assertRaises(Exception, contar_bases,"Affhth_reeet21?*") #testa para seq. invalida
        self.assertEqual(contar_bases("AgtaCGTcAG"), {"A":3, "C":2, "G":3, "T":2}) #maiusculas e minusculas
        self.assertEqual(contar_bases("agtaaCGTcgaG"), {"A":4, "C":2, "G":4, "T":2}) #maiusculas e minusculas
        
    def test_inverso(self):
        self.assertEqual(complemento_inverso("GTACTGATGCAAAACCGGTT"), "AACCGGTTTTGCATCAGTAC") #testa seq normal
        self.assertEqual(complemento_inverso("TCgcGCGCaTgggCC"), "GGCCCATGCGCGCGA") #maiusculas e minusculas
        self.assertEqual(complemento_inverso("tgctgtagtcagtcgtagtcata"), "TATGACTACGACTGACTACAGCA") #minusculas
        self.assertRaises(Exception, complemento_inverso, "aUAaCauaccGg") #testa para RNA
        self.assertRaises(Exception, complemento_inverso, "AuGdUUAa_?/5Caa") #testa para seq. invalida
        
    def test_transcricao(self):
        self.assertEqual(transcricao("TCTATCTCCGGTCCAAGTCA"), "UCUAUCUCCGGUCCAAGUCA")
        self.assertEqual(transcricao("atcgatgatcgtagtctacggtcat"), "AUCGAUGAUCGUAGUCUACGGUCAU")
        self.assertEqual(transcricao("TcTaTCTCCggTCCAagTcA"), "UCUAUCUCCGGUCCAAGUCA")
        self.assertRaises(Exception, transcricao, "aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa") #testa para seq. invalida
        self.assertRaises(Exception, transcricao, "aUAaCauaccGg") #testa para RNA
        
    def test_traducao(self):
        self.assertEqual(traducao("atgtagctgatcgtagctagtctatc"), "M_LIVASL")
        self.assertEqual(traducao("GCAATACCGTGGTCTAGCCTCAAATTCGGCTCTGTTACCTCGAGCGTTATGTGTCAAATGGCG"), "AIPWSSLKFGSVTSSVMCQMA")
        self.assertEqual(traducao("TcTaTCTCCAtgggTCCAagTcAttaA"), "SISMGPSH_") #maiusculas e minusculas
        self.assertRaises(Exception, traducao, "aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa") #testa para seq. invalida
        self.assertRaises(Exception, traducao, "UauaaCuaGucgucuaaUACuccacucGg") #testa para RNA
        
    def test_reading_frames(self):
        self.assertCountEqual(reading_frames("AGTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA"), ['AGTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'TGCGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT', 'GTGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'GCGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT', 'TGATTTTGAAAGGTTTGCGGGGCACAGCTATGACTTGCTTAGCTACGTGTGAGGGAAGGAACTTTTGCGTATTTGTATGTTCACCCGTCTACTACGCA', 'CGTAGTAGACGGGTGAACATACAAATACGCAAAAGTTCCTTCCCTCACACGTAGCTAAGCAAGTCATAGCTGTGCCCCGCAAACCTTTCAAAATCACT'])
        self.assertCountEqual(reading_frames("tgcgggcagattatgtaggttgagagatgcgggagaagttctcgaccttcccgtgggacgtgaacctatcccctaatagagcattccgttcgagcatggc"), ['TGCGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'GCCATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA', 'GCGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'CCATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA', 'CGGGCAGATTATGTAGGTTGAGAGATGCGGGAGAAGTTCTCGACCTTCCCGTGGGACGTGAACCTATCCCCTAATAGAGCATTCCGTTCGAGCATGGC', 'CATGCTCGAACGGAATGCTCTATTAGGGGATAGGTTCACGTCCCACGGGAAGGTCGAGAACTTCTCCCGCATCTCTCAACCTACATAATCTGCCCGCA'])
        self.assertRaises(Exception, reading_frames, "rueth8347y32894u59806y4957813r837gh25gh") #testa para seq. invalida

    def test_get_proteins(self):
        self.assertEqual(get_proteins("ATGACCGTAA"), [])
        self.assertEqual(get_proteins("GCAATACCGTGGTCTAGCCTCAAATTCGGCTCTGTTACCTCGAGCGTTATGTGTCAAATGGCG"), [])
        self.assertEqual(get_proteins("TcTaTCTCCAtgggTCCAagTcAttaA"), ["MGPSH_"])
        self.assertRaises(Exception, get_proteins, "aagtgtgtgtgcgtagctagtctgcattgctagtccatgctatA+76u+*GdUUAa_?/5Caa") #testa para seq. invalida
        self.assertRaises(Exception, get_proteins, "UauaaCuaGucgucuaaUACuccacucGg") #testa para RNA

# TODO: proteinas maiores dois sentidos, proteinas sem_

if __name__ == "__main__":
    unittest.main()