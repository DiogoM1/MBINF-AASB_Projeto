from sequencinator.sequences import reading_frames, traducao, valida
from sequencinator.data.aln_matrices import blossum62
import re


def get_proteins(seq):
    """
    Função que devolve a lista de todas as proteínas ordenadas por tamanho e por ordem alfabética para as do mesmo tamanho

    PARAMETERS
    :param seq: str
    :return:int
    """
    if valida(seq):
        proteinas = []
        for seq in reading_frames(seq):
            proteinas += re.findall(r'(M[A-Z]*_)', traducao(seq))
        sorted_proteinas = sorted(sorted(set(proteinas)), key=len, reverse=True)
        return sorted_proteinas
    else:
        raise Exception('Não é uma sequência de DNA')

# Alignment tools


def replacement_score(letter_a, letter_b):
    """
    Devolve o score do alinhamento do aa A com aa B usando os scores blossum62

    PARAMETERS
    :param letter_a:str
    :param letter_b:str
    :return:int
    """
    return blossum62[letter_a][letter_b]  # Procura na na matriz blossum62 a score, cada letra funciona como coordenada


def seq_aln_score(seq_a, seq_b, spc_cost):
    """
    Devolve o score do alinhamento direto de duas sequencias, com um calculo de scores para espaços uniforme.

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:int
    """
    score = 0 # Inicia o score com valor 0
    for letter in range(0, len(seq_a)):  # Para cada caracter das duas sequências:
        if seq_a[letter] == "-" or seq_b[letter] == "-":  # Se for um espaço dar adicionar o valor que utilizador inseriu para o score de cada espaço
            score += spc_cost
        else:  # Caso contrário quer dizer que ambos os caracteres são letra e como tal usar a blossum62
            score += replacement_score(seq_a[letter].upper(), seq_b[letter].upper())
    return score

