import re


def ler_seq(FileHandle):
    '''
    s
    '''
    seq = ""
    with open(FileHandle, 'r') as file:
        for line in file:
            seq.append(line)
    return seq

# usar regex para validar ????

def ler_FASTA_seq(FileHandle):
    '''
    s
    '''
    seq = ""
    with open(FileHandle, 'r') as file:
        for line in file:

            seq.append(line)
    return seq
