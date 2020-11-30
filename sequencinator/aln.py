from sequencinator.data.aln_matrices import blossum62
from sequencinator.protein import prep_seq_for_aln


def replacement_score(letter_a, letter_b):
    """
    Devolve o score do alinhamento do aa A com aa B usando os scores blossum62

    PARAMETERS
    :param letter_a:str
    :param letter_b:str
    :return:int
    """
    return blossum62[letter_a][letter_b]  # Procura na na matriz blossum62 a score, cada letra funciona como coordenada


def aln_score(seq_a, seq_b, spc_cost):
    """
    Devolve o score do alinhamento direto de duas sequencias, com um calculo de scores para espaços uniforme.

    PARAMETERS
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


def aln_replacement_score_matrix(seq_a, seq_b, spc_cost):
    """
    Devolve o score do alinhamento do aa A com aa B para todas as combinações da matriz usando os scores blossum62

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:matrix
    """
    seq_a, seq_b = prep_seq_for_aln(seq_a), prep_seq_for_aln(seq_b)  # prepara as duas seq
    matrix = []  # inicia as matrizes
    for a in range(0, len(seq_a)):  # inicia todas as linhas da matriz (uma para cada letra da seq A)
        line = []
        for b in range(0, len(seq_b)): # inicia todas as colunas da matriz (uma para cada letra da seq B)
            line.append(aln_score(seq_b[b], seq_a[a], spc_cost))  # para cada posição da matriz descobre o score da substituição dos aa das duas seq
        matrix.append(line)
    return matrix


def aln_origin_traceback(matrix, seq_a, seq_b, start_a, start_b, end_char):
    str_seq_a = ""
    str_seq_b = ""
    a, b = start_a, start_b  # Capaz de aceitar sequências de tamanhos diferentes

    while (a, b) != (0, 0):  # não chegar ao primeiro quadrado e não chegar a um R (reinicio)
        origin = str(matrix[a][b])[0] # Como um alinhamento pode ter várias origens então vai buscar a primeira.
        if origin == end_char:  #e se não não chegar a um R (reinicio)
            break
        elif origin == "D":
            letter_a, letter_b = seq_a[a-1], seq_b[b-1]
            # A e B ficam como as cordenadas no ponto da origem para onde a proxima operação vai ser executada
            a -= 1
            b -= 1
        elif origin == "C":
            letter_a, letter_b = seq_a[a-1], "-"  # Como foi para cima quer dizer que não usou a letra e como tal espaço.
            a -= 1
        elif origin == "E":
            letter_a, letter_b = "-", seq_b[b-1] # Como foi para a esquerda quer dizer que não usou a letra e como tal espaço.
            b -= 1

        # Adição dos aa alinhados para aquele ponto
        str_seq_a += letter_a
        str_seq_b += letter_b
    return [str_seq_a[::-1], str_seq_b[::-1]] # Devolve as duas se invertidas
