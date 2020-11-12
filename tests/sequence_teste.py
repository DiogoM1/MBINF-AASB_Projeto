import unittest

from sequencinator.sequences import valida, contar_bases


class teste_valida_contar_bases(unittest.TestCase):
    
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
        self.assertFalse(valida("AuGdUUAa_?/5Caa"))    
  
    def test_bases(self):
       self.assertEqual(contar_bases("ACGTAG"), {"A":2, "C":1, "G":2, "T":1})
        
    def test_minusculas_b(self):   
       self.assertEqual(contar_bases("atgtttcagaaag"), {"A":5, "C":1, "G":3, "T":4})
        
    def test_nulo(self):        
        self.assertIsNone(contar_bases("AUGGUUUCA"))
        self.assertIsNone(contar_bases("Affhth_reeet21?*"))
        
    def test_maimin(self):
        self.assertEqual(contar_bases("AgtaCGTcAG"), {"A":3, "C":2, "G":3, "T":2})
        self.assertEqual(contar_bases("agtaaCGTcgaG"), {"A":4, "C":2, "G":4, "T":2})
        
        
        
        
if __name__ == "__main__":
    unittest.main()