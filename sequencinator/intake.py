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
    Parser de um ficheiro no formato FASTA para uma String
    Dado um ficheiro aberto, com uma posição de cursor, lê a próxima sequência escrita no ficheiro.
    Retorna uma String com a sequência, ou vazia quando atingir o EOF.
    '''
    seq = ""
    maiores = 0
    while(maiores<2):
        cursor=FileHandle.tell()
        linha=FileHandle.readline()
        if(linha==""): break
        if linha[0]!='>':
            seq+=linha
        else:
            maiores+=1
    FileHandle.seek(cursor)
    return seq.replace("\n", "")

# main

file = open(r"D:\Universidade\Mestradooo\Algoritmos para Análise de Sequências Biológicas\MBINF-AASB_Projeto\sequencinator\protein.txt", mode='r', encoding='utf-8')
todas_as_seq = []
i=0
while(True):
    i+=1
    input("Carregue no ENTER para ler uma sequência (interromper com Ctrl+C)")
    seq=ler_FASTA_seq(file)
    if seq!="":
        todas_as_seq.append(seq)
        print("Sequencia", i)
        print(seq)
    else: 
        break
file.close

