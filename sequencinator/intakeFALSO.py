def ler_seq2(FileHandle):
    '''
    Função que recebe um ficheiro aberto e que devolve uma sequência (DNA, RNA e aminoácidos)
    PARAMETERS
    seq: str
        a sequência de DNA, RNA e aminoácidos
    returns: sequência contida no ficheiro
    '''
    seq=""; char='';TAMANHO_MINIMO=5
    while(True):
        char=FileHandle.read(1).upper()
        if (char=='A' or char=='T' or char=='C' or char=='G'): seq+=char
        elif (char=='\n'): # para sequências em múltiplas linhas
            i=0 # empata-nabos
        elif (char==''): # para o EOF
            if len(seq)>=TAMANHO_MINIMO: return seq
            else: break
        else:
            if len(seq)>=TAMANHO_MINIMO: return seq
            else: seq=""
            #senão é pequeno demais, continua a procurar


def ler_seq(FileHandle):
    '''
    Função que recebe um ficheiro aberto e que devolve uma sequência (DNA, RNA e aminoácidos)
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
    PARAMETERS
    seq: str
        a sequência de DNA, RNA e aminoácidos
    returns: sequência contida no ficheiro FASTA
    '''
    seq = ""
    maiores = 0 # quantas vezes apareceu o símbolo '>', usado para travar o cíclo no fim da sequência
    while(maiores<2):
        cursor=FileHandle.tell()
        linha=FileHandle.readline().upper()
        if(linha==""): break
        if linha[0]!='>':
            seq+=linha
        else: # if linha[0]=='>'
            maiores+=1
    FileHandle.seek(cursor)
    return seq.replace("\n", "")

###### falta conseguir diferenciar seq de DNA, RNA e de aminoacidos. matem-me

################
##### MAIN #####
################

todas_as_seq = []; i=0; seq=""
path = input("indique o path do ficheiro que quer abrir:" )
file = open(path, mode='r', encoding='utf-8')
# pequeno input/menu para escolher que função chamar
# isto está marretado e podia-se serializar melhor o código para evitar correr tanto if, safoda :shrug:
print("1 - Não FASTA\n2 - FASTA")
while(True):
    op=input()
    if op=="1" or op=="2": break 

while(True):
    i+=1
    input("Carregue no ENTER para ler uma sequência (interromper com Ctrl+C)")
    if(op=="1"): seq=ler_seq(file)
    if(op=="2"): seq=ler_FASTA_seq(file)
    if seq!="":
        todas_as_seq.append(seq)
        print("Sequencia", i)
        print(seq)
    else: 
        print(seq)
        break
file.close