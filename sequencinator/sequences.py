from .data.dic_codao2amino import gene_code as gene_code
import re

def complemento_inverso(seq):
    '''
    '''
    return seq[::-1].upper().replace('T', 'a').replace('A', 't').replace('G', 'c').replace('C', 'g').upper()


def transcricao(seq):
    '''
    '''
    return complemento_inverso(transcricao(seq)).replace('T', 'U')


def traducao(seq):
    '''
    '''
    amino = ''
    for i in range(0, len(seq), 3):
        codao = seq.upper()[i: i + 3]
        if len(codao) == 3:
            amino += gene_code[codao]
    return amino


def reading_frames(seq):
    '''
    '''
    orf = []
    for i in range(0, 3):
        orf.append[seq.upper()[i:], complemento_inverso(seq.upper())[i:]]
    return orf


def valida(seq):
    '''
    Verifica se a sequência de ADN é válida
   
    Parameters
        ----------
        seq : str
            Sequência de ADN  

   Returns
        ------
        True ou False
 
    '''
    seq = seq.upper()     
    if  re.search("[^TAGC]",seq) == None:        
        validade = True
    else:
        validade = False
    return validade

