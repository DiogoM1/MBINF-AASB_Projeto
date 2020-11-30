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


def prep_seq_for_aln(seq):
    """
    Função de limpeza e formatação para as seq que vão ser alinhas por matriz

    PARAMETERS
    :param seq:str
    :return:str
    """
    seq = seq.strip().upper() # Remove espaços " " antes e depois e garante que a seq está em maiusculas
    if seq[0] != "-":  # adiciona o espaço, representado por "-", necessário à formação de matrizes de alinhamento, no inicio
        seq = "-" + seq
    return seq


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


def aln_needleman(seq_a, seq_b, spc_cost):
    """
    Devolve a matriz de alinhamento do algoritmo Needleman-Wunsch

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:matrix
    """
    matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)
    for a in range(0, len(matrix)):
        for b in range(0, len(matrix[0])):
            if a == 0 and b == 0:  # ponto de cordenadas (0,0)
                merito = 0
            elif a == 0:  # toda a primeira linha
                merito = matrix[a][b-1]+spc_cost
            elif b == 0:  # toda a primeira coluna
                merito = matrix[a-1][b]+spc_cost
            else:  # restantes campos da matriz
                diag = matrix[a-1][b-1]+matrix[a][b]  # soma do merito de alinhamento de origem diagonal
                vert = matrix[a-1][b]+spc_cost  # soma do merito de alinhamento de origem vertical
                hori = matrix[a][b-1]+spc_cost  # soma do merito de alinhamento de origem horizontal
                merito = max(diag, vert, hori)  # seleção da score do melhor alinhamento possivel
            matrix[a][b] = merito
    return matrix


def aln_needleman_best(seq_a, seq_b, spc_cost):
    """
    Devolve o score do melhor alinhamento de duas sequencias com o algoritmo Needleman-Wunsch

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:int
    """
    return aln_needleman(seq_a, seq_b, spc_cost)[-1][-1]  # seleção do ultimo campo da matriz, correspondente ao score total do melhor alinhamento


def aln_needleman_origin(seq_a, seq_b, spc_cost):
    """
    Devolve a matriz da origem do melhor alinhamento do algoritmo Needleman-Wunsch

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:
    """
    rep_matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)  # rep = replacement
    needleman_matrix = aln_needleman(seq_a, seq_b, spc_cost)
    origin_matrix = [[0 for b in range(0, len(needleman_matrix[0]))] for a in range(0, len(needleman_matrix))] # iniciar uma nova matriz

    for a in range(0, len(needleman_matrix)):
        for b in range(0, len(needleman_matrix[0])):
            origin = ""
            if a == 0 and b == 0:
                origin += ""
            elif a == 0:
                origin += "E"
            elif b == 0:
                origin += "C"
            else:
                cell_max = needleman_matrix[a][b]  # consultar na matriz de needleman o valor do melhor alinhamento para aquele ponto
                # Comparação com as 3 origens possiveis para perceber a origem do valor da matriz de needleman
                if cell_max == needleman_matrix[a-1][b-1]+rep_matrix[a][b]:
                    origin += "D"  # Diagonal
                if cell_max == needleman_matrix[a][b-1]+int(spc_cost):
                    origin += "E"  # Esquerda/horizontal
                if cell_max == needleman_matrix[a-1][b]+int(spc_cost):
                    origin += "C"  # Cima/vertical
            origin_matrix[a][b] = origin  # adição das letras que codificam a orogem à nova matriz
    return origin_matrix


def aln_needleman_traceback(seq_a, seq_b, spc_cost):
    """
    Devolve as duas seq no seu melhor alinhamento com o algoritmo de needleman-wunsch

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return: list[str, str]
    """
    matrix = aln_needleman_origin(seq_a, seq_b, spc_cost)  # Utiliza a matriz de origens e corre esta do fim para o inicio de forma dinámica
    str_seq_a = ""
    str_seq_b = ""
    a, b = len(seq_a), len(seq_b)  # Capaz de aceitar sequências de tamanhos diferentes

    while (a, b) != (0, 0):  # não chegar ao primeiro quadrado
        origin = str(matrix[a][b])[0] # Como um alinhamento pode ter várias origens então vai buscar a primeira.
        if origin == "D":
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