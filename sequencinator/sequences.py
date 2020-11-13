from .data.dic_codao2amino import gene_code as gene_code
import re


def valida_rna(seq):
    '''
    Verifica se a sequência de RNA é válida.

    Parameters
        ----------
        seq : str
            Sequência de ADN.

    Returns
        ------
        True ou False

    '''
    seq = seq.upper()
    if re.search(r'[^UAGC]', seq) is None:
        validade = True
    else:
        validade = False
    return validade


def valida(seq):
    '''
    Verifica se a sequência de ADN é válida.

    Parameters
        ----------
        seq : str
            Sequência de ADN.

    Returns
        ------
        True ou False

    '''
    seq = seq.upper()
    if re.search(r'[^TAGC]', seq) is None:
        validade = True
    else:
        validade = False
    return validade


def contar_bases(seq):
    '''
    Conta as bases de uma sequência de ADN.

    Parameters
    ----------
    seq : str
        Sequência de ADN.

    Returns
    -------
    Dicionário com a contagem.

    '''
    if valida(seq) or valida_rna(seq):
        pb = {}
        seq = seq.upper()
        for base in seq:
            pb[base] = pb.get(base, 0) + 1
        return pb
    else:
        raise Exception('Não é uma sequência de DNA')


def complemento_inverso(seq):
    '''
    Função que devolve o complemento inverso de uma sequência de DNA
    PARAMETERS
    seq: str
        a sequência de
    returns: str do complemento inverso da sequencia fornecida
    '''
    if valida(seq):
        return seq[::-1].upper().replace('T', 'a').replace('A', 't').replace('G', 'c').replace('C', 'g').upper()
    else:
        raise Exception('Não é uma sequência de DNA')


def transcricao(seq):
    '''
    Função que devolve a transcrição de uma sequência de DNA

    PARAMETERS
    seq: str
        a sequência de
    returns: str do mRNA para a sequencia fornecida
    '''
    if valida(seq):
        return complemento_inverso(complemento_inverso(seq)).replace('T', 'U')
    else:
        raise Exception('Não é uma sequência de DNA')


def traducao(seq):
    '''
    Função que devolve a tradução de uma sequência de DNA

    PARAMETERS
    seq: str
        a sequência de
    returns: str de cadeia de aminoacidos
    '''
    if valida(seq):
        amino = ''
        for i in range(0, len(seq), 3):
            codao = seq.upper()[i: i + 3]
            if len(codao) == 3:
                amino += gene_code[codao]
        return amino
    else:
        raise Exception('Não é uma sequência de DNA')


def reading_frames(seq):
    '''
    Função que devolve uma lista com as reading frames

    PARAMETERS
    seq: str
        a sequência de
    returns: list das diferentes ORFs
    '''
    if valida(seq):
        orf = []
        for i in range(0, 3):
            orf.extend((seq.upper()[i:], complemento_inverso(seq.upper())[i:]))
        return orf
    else:
        raise Exception('Não é uma sequência de DNA')
