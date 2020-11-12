def ler_seq(FileHandle):
    '''
    Função que recebe um ficheiro aberto e que devolve uma sequência (DNA, RNA e aminoácidos)
    ------------------------------------
    PARAMETERS
    seq: str
        a sequência de DNA, RNA e aminoácidos
    returns: sequência contida no ficheiro
    '''
    seq=""
    while(True):
        char=FileHandle.read(1).upper()
        if (char>='A' and char <='Z'): seq+=char
        if char=="": return seq



def ler_FASTA_seq(FileHandle):
    '''
    Função que recebe um ficheiro FASTA aberto e que devolve uma sequência (DNA, RNA e aminoácidos)
    ------------------------------------
    PARAMETERS
    seq: str
        a sequência de DNA, RNA e aminoácidos
    returns: sequência contida no ficheiro FASTA
    '''
    seq = ""
    maiores = 0
    while(maiores<2):
        cursor=FileHandle.tell()
        linha=FileHandle.readline().upper()
        if(linha==""): break
        if linha[0]!='>':
            seq+=linha
        else:
            maiores+=1
    FileHandle.seek(cursor)
    return seq.replace("\n", "")
