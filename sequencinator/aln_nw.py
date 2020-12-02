from sequencinator import aln_replacement_score_matrix
from sequencinator.aln import aln_origin_traceback


def aln_pd_nw(matrix, spc_cost):
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


def aln_nw(seq_a, seq_b, spc_cost):
    """
    Devolve a matriz de alinhamento do algoritmo Needleman-Wunsch

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:matrix
    """
    base_matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)
    final_matrix = aln_pd_nw(base_matrix, spc_cost)
    return final_matrix


def aln_nw_best(seq_a, seq_b, spc_cost):
    """
    Devolve o score do melhor alinhamento de duas sequencias com o algoritmo Needleman-Wunsch

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:int
    """
    return aln_nw(seq_a, seq_b, spc_cost)[-1][-1]  # seleção do ultimo campo da matriz, correspondente ao score total do melhor alinhamento


def aln_nw_origin(seq_a, seq_b, spc_cost):
    """
    Devolve a matriz da origem do melhor alinhamento do algoritmo Needleman-Wunsch

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:
    """
    rep_matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)  # rep = replacement
    needleman_matrix = aln_nw(seq_a, seq_b, spc_cost)
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
                if cell_max == needleman_matrix[a-1][b]+int(spc_cost):
                    origin += "C"  # Cima/vertical
                if cell_max == needleman_matrix[a][b-1]+int(spc_cost):
                    origin += "E"  # Esquerda/horizontal
            origin_matrix[a][b] = origin  # adição das letras que codificam a orogem à nova matriz
    return origin_matrix


def aln_nw_traceback(seq_a, seq_b, spc_cost):
    """
    Devolve as duas seq no seu melhor alinhamento com o algoritmo de needleman-wunsch

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return: list[str, str]
    """
    matrix = aln_nw(seq_a, seq_b, spc_cost)
    origin_matrix = aln_nw_origin(seq_a, seq_b, spc_cost)
    aln = aln_origin_traceback(origin_matrix, seq_a, seq_b, len(seq_a), len(seq_b), "")  # Começa no fim da matrix
    return aln