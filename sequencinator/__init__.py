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
aln_replacement_score_matrix - Devolve o score do alinhamento do aa A com aa B para todas as combinações da matriz usando os scores blossum62

Alignment
replacement_score - Devolve o score do alinhamento do aa A com aa B usando os scores blossum62
aln_score - Devolve o score do alinhamento direto de duas sequencias, com um calculo de scores para espaços uniforme.

    Needleman-Wunsch Algorithm
    aln_nw - Devolve a matriz de alinhamento do algoritmo Needleman-Wunsch
    aln_nw_origin - Devolve o score do melhor alinhamento de duas sequencias com o algoritmo Needleman-Wunsch
    aln_nw_traceback - Devolve as duas seq no seu melhor alinhamento com o algoritmo de Needleman-Wunsch

    Smith-Waterman Algorithm
    aln_sw - Devolve a matriz de alinhamento do algoritmo Smith-Waterman
    aln_sw_origin - Devolve o score do melhor alinhamento de duas sequencias com o algoritmo Smith-Waterman
    aln_sw_traceback - Devolve as duas seq no seu melhor alinhamento com o algoritmo de Smith-Waterman

"""

__version__ = "1.0"
__author__ = "André Ramos, Angelina Eiras, Diogo Macedo, Pedro Gonçalves"
__license__ = ""

from .intake import ler_FASTA_seq, ler_seq
from .sequences import complemento_inverso, transcricao, traducao, reading_frames, valida_rna, valida, contar_bases
from .protein import get_proteins
from .aln import replacement_score, aln_score, aln_replacement_score_matrix
from .aln_nw import aln_nw, aln_nw_origin, aln_nw_traceback
from .aln_sw import aln_sw, aln_sw_origin, aln_sw_traceback
