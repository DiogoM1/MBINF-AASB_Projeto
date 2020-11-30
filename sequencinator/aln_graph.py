def aln_graph(seq_a, seq_b):
    # http://www.bioinformatics.uni-muenster.de/teaching/courses-2016/bioinf1/download/bioinfo1_3_2016_handouts.pdf
    # http://profs.scienze.univr.it/~liptak/FundBA/slides/PWAl_inPractice.pdf
    # http://www.bpc.uni-frankfurt.de/guentert/wiki/images/5/58/121102_PairwiseAlignment.pdf
    dot_matrix = []
    seq_a, seq_b = seq_a.upper(), seq_b.upper()
    for a in range(len(seq_a)):
        dot_matrix.append("")
        for b in range(len(seq_b)):
            if seq_a[a] == seq_b[b]:
                dot_matrix[a] += "*"
            else:
                dot_matrix[a] += "_"
    return dot_matrix