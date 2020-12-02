from sequencinator import aln_replacement_score_matrix
from sequencinator.aln import aln_origin_traceback
from sequencinator.matrix_tools import find_last_max, find_all_max


def aln_pd_sw(matrix, spc_cost):
    for a in range(0, len(matrix)):
        for b in range(0, len(matrix[0])):
            if a == 0 and b == 0:  # ponto de cordenadas (0,0)
                merito = 0
            elif a == 0:  # toda a primeira linha
                merito = max(matrix[a][b-1]+spc_cost,0)
            elif b == 0:  # toda a primeira coluna
                merito = max(matrix[a-1][b]+spc_cost,0)
            else:  # restantes campos da matriz
                diag = matrix[a-1][b-1]+matrix[a][b]  # soma do merito de alinhamento de origem diagonal
                vert = matrix[a-1][b]+spc_cost  # soma do merito de alinhamento de origem vertical
                hori = matrix[a][b-1]+spc_cost  # soma do merito de alinhamento de origem horizontal
                merito = max(diag, vert, hori, 0)  # seleção da score do melhor alinhamento possivel ou 0 (reinicia)
            matrix[a][b] = merito
    return matrix


def aln_sw(seq_a, seq_b, spc_cost):
    """
    Devolve a matriz de alinhamento do algoritmo Smith-Waterman

    PARAMETERS
    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return:matrix
    """
    base_matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)
    final_matrix = aln_pd_sw(base_matrix, spc_cost)
    return final_matrix


def aln_sw_origin(seq_a, seq_b, spc_cost):
    rep_matrix = aln_replacement_score_matrix(seq_a, seq_b, spc_cost)  # rep = replacement
    sw_matrix = aln_sw(seq_a, seq_b, spc_cost)
    origin_matrix = [[0 for b in range(0, len(sw_matrix[0]))] for a in range(0, len(sw_matrix))] # iniciar uma nova matriz

    for a in range(0, len(sw_matrix)):
        for b in range(0, len(sw_matrix[0])):
            origin = ""
            if a == 0 or b == 0:
                origin += "R"
            else:
                cell_max = sw_matrix[a][b]  # consultar na matriz de needleman o valor do melhor alinhamento para aquele ponto
                # Comparação com as 3 origens possiveis para perceber a origem do valor da matriz de needleman
                if cell_max == 0:
                    origin += "R"
                else:
                    if cell_max == sw_matrix[a-1][b-1]+rep_matrix[a][b]:
                        origin += "D"  # Diagonal
                    if cell_max == sw_matrix[a][b-1]+int(spc_cost):
                        origin += "E"  # Esquerda/horizontal
                    if cell_max == sw_matrix[a-1][b]+int(spc_cost):
                        origin += "C"  # Cima/vertical
            origin_matrix[a][b] = origin  # adição das letras que codificam a orogem à nova matriz
    return origin_matrix


def aln_sw_traceback(seq_a, seq_b, spc_cost):
    """
    Devolve as duas seq no seu melhor alinhamento com o algoritmo de smith-waterman

    :param seq_a:str
    :param seq_b:str
    :param spc_cost:int (can be negative)
    :return: list[str, str]
    """
    matrix = aln_sw(seq_a, seq_b, spc_cost)
    max_indexes = find_all_max(matrix)
    origin_matrix = aln_sw_origin(seq_a, seq_b, spc_cost)
    aln = []
    for origin in max_indexes:
        a, b = origin[0], origin[1]
        aln.append(aln_origin_traceback(origin_matrix, seq_a, seq_b, a, b, "R"))  # Começa no maior ultimo valor da matrix
    return aln
