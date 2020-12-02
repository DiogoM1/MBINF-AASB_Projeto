def max_matrix(matrix):
    # finds max value
    # https://www.educba.com/numpy-argmax/
    max_list = []
    for a in range(0, len(matrix)):
        max_list.append(max(matrix[a]))
    return max(max_list)


def find_last_max(matrix):
    # finds the last max coordinates
    max = max_matrix(matrix)
    max_index_sum = 0
    index_a, index_b = 0, 0
    for a in range(0, len(matrix)):
        for b in range(0, len(matrix[0])):
            # TODO: fix "desempate" logic
            if matrix[a][b] == max and a+b >= max_index_sum:
                index_a, index_b = a, b
    return index_a, index_b


def find_all_max(matrix):
    # finds the all max coordinates
    max = max_matrix(matrix)
    index_list = []
    for a in range(0, len(matrix)):
        for b in range(0, len(matrix[0])):
            if matrix[a][b] == max:
                index_list.append((a, b))
    return index_list