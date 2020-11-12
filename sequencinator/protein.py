from .sequences import reading_frames, traducao
import re


def get_proteins(seq):
    '''
    Função que devolve a lista de todas as proteínas ordenadas por tamanho e por ordem alfabética para as do mesmo tamanho

    PARAMETERS
    seq: str
        a sequência de
    returns: lista de proteinas
    '''
    orf = reading_frames(seq)
    proteinas = []
    for seq in orf:
        proteinas += re.findall(r'(M[A-Z]*_)', traducao(seq))
    sorted_proteinas = sorted(sorted(set(proteinas)), key=len, reverse=True)
    return sorted_proteinas
