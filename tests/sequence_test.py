import unittest

from sequencinator.sequences import complemento_inverso, transcricao, traducao, reading_frames, valida, contar_bases, \
    valida_rna
from sequencinator.protein import get_proteins, prep_seq_for_aln
from sequencinator.aln_nw import aln_nw_best
from sequencinator.aln_sw import aln_sw, aln_sw_origin, aln_sw_traceback
from sequencinator import replacement_score, aln_score, aln_replacement_score_matrix, aln_nw, \
    aln_nw_origin, aln_nw_traceback
from sequencinator.matrix_tools import max_matrix, find_last_max


class TesteFuncoes(unittest.TestCase):
    
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
        self.assertRaises(Exception, contar_bases, "Affhth_reeet21?*") #testa para seq. invalida
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

    def test_replacement_score(self):
        self.assertEqual(replacement_score("A", "A"), 4)

# TODO: Melhores testes para todas as da semana 5

    def test_seq_aln_score(self):
        self.assertEqual(aln_score("-HGWAG", "PHSW-G", -8), 9)

    def test_prep_seq_for_al(self):
        self.assertEqual(prep_seq_for_aln("-HGWAG"), "-HGWAG")
        self.assertEqual(prep_seq_for_aln("HGWAG"), "-HGWAG")

    def test_aln_replacement_score_matrix(self):
        self.assertEqual(aln_replacement_score_matrix(" PHSWG", "HGWAG", -8), [[-8, -8, -8, -8, -8, -8],
                                                                  [-8, -2, -2, -4, -1, -2],
                                                                  [-8, 8, -2, -2, -2, -2],
                                                                  [-8, -1, 0, -3, 1, 0],
                                                                  [-8, -2, -2, 11, -3, -2],
                                                                  [-8, -2, 6, -2, 0, 6]])
        self.assertEqual(aln_replacement_score_matrix("PHSWG", "HGWAG", -8), [[-8, -8, -8, -8, -8, -8],
                                                                  [-8, -2, -2, -4, -1, -2],
                                                                  [-8, 8, -2, -2, -2, -2],
                                                                  [-8, -1, 0, -3, 1, 0],
                                                                  [-8, -2, -2, 11, -3, -2],
                                                                  [-8, -2, 6, -2, 0, 6]])
        self.assertEqual(aln_replacement_score_matrix("PHSWG", "HGwAG", -8), [[-8, -8, -8, -8, -8, -8],
                                                                  [-8, -2, -2, -4, -1, -2],
                                                                  [-8, 8, -2, -2, -2, -2],
                                                                  [-8, -1, 0, -3, 1, 0],
                                                                  [-8, -2, -2, 11, -3, -2],
                                                                  [-8, -2, 6, -2, 0, 6]])

    def test_aln_needleman(self):
        """
        """

    def test_aln_needleman_origin(self):
        self.assertEqual(aln_nw_origin('IPGKASYD', 'VSPAGMASGYD', -4), [
            ["", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E"],
            ["C", "D", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E"],
            ["C", "C", "D", "D", "E", "E", "E", "E", "E", "E", "E", "E"],
            ["C", "C", "D", "C", "D", "D", "E", "E", "E", "DE", "E", "E"],
            ["C", "C", "DC", "DC", "C", "DC", "D", "DE", "D", "E", "E", "E"],
            ["C", "C", "D", "DC", "D", "D", "DC", "D", "E", "E", "E", "E"],
            ["C", "C", "D", "D", "C", "D", "D", "C", "D", "E", "E", "E"],
            ["C", "C", "C", "D", "C", "C", "D", "C", "C", "D", "D", "E"],
            ["C", "C", "C", "D", "C", "C", "C", "DC", "C", "D", "C", "D"]
        ])

    def test_aln_needleman_best(self):
        self.assertEqual(aln_nw_best('HGWAG', 'PHSWG', -8), 9)
        self.assertEqual(aln_nw_best(' HGWAG ', ' PHSWG', -8), 9)
        self.assertEqual(aln_nw_best('HGwAG', 'phswg', -8), 9)
        self.assertEqual(aln_nw_best('HGWAG', 'PHSWG', -8), 9)
        self.assertEqual(aln_nw_best('IPGKASYD', 'VSPAGMASGYD', -4), 24)

    def test_aln_needleman_traceback(self):
        self.assertEqual(aln_nw_traceback('HGWAG', 'PHSWG', -8), ["-HGWAG", "PHSW-G"])
        self.assertEqual(aln_nw_traceback('IPGKASYD', 'VSPAGMASGYD', -4), ["I-P-GKAS-YD", "VSPAGMASGYD"])

    def test_find_last_max(self):
        self.assertEqual(find_last_max([[-8, -8, -8, -8, -8, -8],
                                        [-8, -2, -2, -4, -1, -2],
                                        [-8, 8, -2, -2, -2, -2],
                                        [-8, -1, 0, -3, 1, 0],
                                        [-8, -2, -2, 11, -3, -2],
                                        [-8, -2, 6, -2, 0, 6]]), (4, 3))

    def test_max_matrix(self):
        self.assertEqual(max_matrix([[-8, -8, -8, -8, -8, -8],
                                        [-8, -2, -2, -4, -1, -2],
                                        [-8, 8, -2, -2, -2, -2],
                                        [-8, -1, 0, -3, 1, 0],
                                        [-8, -2, -2, 11, -3, -2],
                                        [-8, -2, 6, -2, 0, 6]]), 11)

    def test_aln_sw(self):
        self.assertEqual(aln_sw(" PHSWG", "HGWAG", -8),
                         [[0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0],
                          [0, 8, 0, 0, 0, 0],
                          [0, 0, 8, 0, 1, 0],
                          [0, 0, 0, 19, 11, 3],
                          [0, 0, 6, 11, 19, 17]])

    def test_aln_sw_origin(self):
        self.assertEqual(aln_sw_origin(" PHSWG", "HGWAG", -8),
                         [["R", "R", "R", "R", "R", "R"],
                          ["R", "R", "R", "R", "R", "R"],
                          ["R", "D", "R", "R", "R", "R"],
                          ["R", "R", "D", "R", "D", "R"],
                          ["R", "R", "R", "D", "E", "E"],
                          ["R", "R", "D", "C", "D", "D"]])

    def test_aln_sw_traceback(self):
        self.assertEqual(aln_sw_traceback('HGWAG', 'PHSWG', -8), ["HGWA", "HSWG"])


if __name__ == "__main__":
    unittest.main()