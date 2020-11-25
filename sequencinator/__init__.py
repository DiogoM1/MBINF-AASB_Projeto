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
aln_score - Devolve o score do alinhamento direto de duas sequencias, com um calculo de scores para espaços uniforme.
aln_replacement_score_matrix - Devolve o score do alinhamento do aa A com aa B para todas as combinações da matriz usando os scores blossum62
aln_needleman - Devolve a matriz de alinhamento do algoritmo Needleman-Wunsch
aln_needleman_origin - Devolve o score do melhor alinhamento de duas sequencias com o algoritmo Needleman-Wunsch
aln_needleman_traceback - Devolve as duas seq no seu melhor alinhamento com o algoritmo de needleman-wunsch

"""

__version__ = "1.0"
__author__ = "André Ramos, Angelina Eiras, Diogo Macedo, Pedro Gonçalves"
__license__ = ""

from .sequences import complemento_inverso, transcricao, traducao, reading_frames, valida_rna, valida, contar_bases
from .protein import get_proteins, replacement_score, aln_score, aln_replacement_score_matrix, aln_needleman, aln_needleman_origin, aln_needleman_traceback
from .intake import ler_FASTA_seq, ler_seq
