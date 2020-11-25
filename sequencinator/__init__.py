# -*- coding: utf-8 -*-
"""
Conjunto de módulos que permitem trabalhar com sequencias biológias

Intake
ler_FASTA_seq - Função que recebe um ficheiro FASTA aberto e que devolve uma sequência (DNA, RNA e aminoácidos)
ler_seq - Função que recebe um ficheiro aberto e que devolve uma sequência (DNA, RNA e aminoácidos)

RNA/DNA
valida_rna - Verifica se a sequência de RNA é válida.
valida - Verifica se a sequência de ADN é válida.
contar_bases - Conta as bases de uma sequência de ADN.
complemento_inverso - Função que devolve o complemento inverso de uma sequência de DNA
transcricao - Função que devolve a transcrição de uma sequência de DNA
traducao - Função que devolve a tradução de uma sequência de DNA
reading_frames - Função que devolve uma lista com as reading frames

Proteinas
get_proteins - Função que devolve a lista de todas as proteínas ordenadas por tamanho e por ordem alfabética para as do mesmo tamanho
replacement_score - Devolve o score do alinhamento do aa A com aa B usando os scores blossum62
seq_aln_score - Devolve o score do alinhamento direto de duas sequencias, com um calculo de scores para espaços uniforme.

"""

__version__ = "1.0"
__author__ = "André Ramos, Angelina Eiras, Diogo Macedo, Pedro Gonçalves"
__license__ = ""

from .sequences import complemento_inverso, transcricao, traducao, reading_frames, valida_rna, valida, contar_bases
from .protein import get_proteins, replacement_score, seq_aln_score
from .intake import ler_FASTA_seq, ler_seq
