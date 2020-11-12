def ler_seq(FileHandle):
    '''
    Parser de uma sequência arbitrária de caracteres
    É considerado qualquer conjunto de símbolos ACTG de comprimento superior a TAMANHO_MINIMO,
    que seja contíguo ou separado por '\n' (para apanhar sequências multi-linha) 
    '''
    seq=""; char='';TAMANHO_MINIMO=5
    while(True):
        char=file.read(1).upper()
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

###### falta conseguir diferenciar seq de DNA, RNA e de aminoacidos. matem-me

def ler_FASTA_seq(FileHandle):
    '''
    Parser de um ficheiro no formato FASTA para uma String
    Dado um ficheiro aberto, com uma posição de cursor, lê a próxima sequência escrita no ficheiro.
    Retorna uma String com a sequência, ou vazia quando atingir o EOF.
    '''
    seq = ""
    maiores = 0 # quantas vezes apareceu o símbolo '>', usado para travar o cíclo no fim da sequência
    while(maiores<2):
        cursor=FileHandle.tell()
        linha=FileHandle.readline()
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
path = input("indique o path do ficheiro FASTA que quer abrir:" )
file = open(path, mode='r', encoding='utf-8')
# pequeno input/menu para escolher que função chamar
# isto está marretado e podia-se serializar melhor o código para evitar correr tanto if, safoda :shrug:
print("1 - FASTA\n2 - Não FASTA")
while(True):
    op=input()
    if op=="1" or op=="2": break 

while(True):
    i+=1
    input("Carregue no ENTER para ler uma sequência (interromper com Ctrl+C)")
    if(op=="1"): seq=ler_FASTA_seq(file)
    if(op=="2"): seq=ler_seq(file); 
    if seq!="":
        todas_as_seq.append(seq)
        print("Sequencia", i)
        print(seq)
    else: 
        break
file.close