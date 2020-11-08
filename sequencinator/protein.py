from .sequences import reading_frames, traducao
import re


def reading_frames(seq):
    '''
    '''
    orf = reading_frames(seq)
    proteinas = []
    for seq in orf:
        proteinas += re.findall(r'(M[A-Z]*_)', traducao(seq))
    sorted_proteinas = sorted(sorted(set(proteinas)), key=len, reverse=True)
    return sorted_proteinas
