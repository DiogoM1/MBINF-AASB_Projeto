from sequencinator.sequences import reading_frames, traducao, valida
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


