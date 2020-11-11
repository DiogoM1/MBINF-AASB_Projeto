from .data.dic_codao2amino import gene_code as gene_code


def complemento_inverso(seq):
    '''
    Função que devolve o complemento inverso de uma sequência de DNA
    PARAMETERS
    seq: str
        a sequência de
    returns: str do complemento inverso da sequencia fornecida
    '''
    return seq[::-1].upper().replace('T', 'a').replace('A', 't').replace('G', 'c').replace('C', 'g').upper()


def transcricao(seq):
    '''
    Função que devolve a transcrição de uma sequência de DNA

    PARAMETERS
    seq: str
        a sequência de
    returns: str do mRNA para a sequencia fornecida
    '''
    return complemento_inverso(transcricao(seq)).replace('T', 'U')


def traducao(seq):
    '''
    Função que devolve a tradução de uma sequência de DNA

    PARAMETERS
    seq: str
        a sequência de
    returns: str de cadeia de aminoacidos
    '''
    amino = ''
    for i in range(0, len(seq), 3):
        codao = seq.upper()[i: i + 3]
        if len(codao) == 3:
            amino += gene_code[codao]
    return amino


def reading_frames(seq):
    '''
    Função que devolve uma lista com as reading frames

    PARAMETERS
    seq: str
        a sequência de
    returns: lista das diferentes ORFs
    '''
    orf = []
    for i in range(0, 3):
        orf.append[seq.upper()[i:], complemento_inverso(seq.upper())[i:]]
    return orf
