"""Swap and euler algorithms.

Shuffling biological sequences
- https://www.sciencedirect.com/science/article/pii/S0166218X97814564
"""

# from constants import S1
import random
from utils import same_klets
from altschul import edge_ordering
from pprint import pprint
import random


def is_k_cyclic(seq: str, k: int) -> bool:
    """Does no check if k makes sense."""
    return seq[: (k - 1)] == seq[-(k - 1) :]


def random_rotation(seq: str, k: int, *, m: int | None = None) -> str:
    """Preserves k-lets.

    Note that this is not a proper rotation. It does not
    conserve the characters of the sequence.

    m can be passed as an int for testing.
    """
    assert is_k_cyclic(seq, k)

    n = len(seq)
    if m is None:
        m = random.randint(k, n)
    assert k <= m <= n

    rotated = seq[m - 1 :] + seq[k - 1 : m]
    for i in range(m + 1, m + k - 1):
        idx = i if i <= n else i % n + k - 1
        rotated += seq[idx - 1]

    assert same_klets(seq, rotated, k)

    return rotated


def back_ran_walk(einv, vertex, T, beta):
    # TODO:better to do once instead of each time the function it is call
    random_element = random.choice(einv[vertex])
    # einv[vertex].remove(random_element)

    if not T[random_element[1:]] and random_element[1:] != beta:
        T[random_element[1:]].append(random_element)

    return random_element[1:]


def random_seq(eord, vertex, T):
    if eord[vertex]:
        # TODO:better to do once instead of each time the function it is call
        random_element = random.choice(eord[vertex])
        eord[vertex].remove(random_element)
    elif T[vertex]:
        random_element = T[vertex][0]
        T[vertex].remove(random_element)
    else:
        raise "Error"

    return random_element[1:]


def victor_algorithm(seq: str, k: int):
    """victor algorithm for generating a uniform random permutation of seq."""
    len_seq = len(seq)
    assert k < len_seq

    # 1. Construct D_k
    eord = edge_ordering(seq, k)
    einv = edge_ordering(seq[::-1], k)

    # 2. If cyclic, random rotation.
    if is_k_cyclic(seq, k):
        rot = random_rotation(seq, k)
        fst = rot[0 : k - 1]
        lst = rot[0 : k - 1]
    else:
        # If acyclic, add dummy edge.
        fst = seq[: k - 1]
        lst = seq[-(k - 1) :]

    # 3.
    inv_T = {i: [] for i in einv.keys()}
    vertex = lst[::-1]
    flag = True
    step = 0
    while flag:  # change into raise error
        # print("vertex:",vertex)
        # print("T:",inv_T,lst[::-1] )
        # print("einv",einv)
        vertex = back_ran_walk(einv, vertex, inv_T, lst[::-1])
        cond = [inv_T[vertex] for vertex in einv.keys() if vertex != lst[::-1]]
        # print("cond",cond)
        flag = not all(cond)
        step += 1
        # if step==len_seq:
        #     break #should be rise error

    # b
    T = {}
    for i in inv_T.keys():
        for value in inv_T[i]:
            edge = value[::-1]
            node = edge[: k - 1]
            if node in T:
                T[node].append(edge)
            else:
                T[node] = [edge]

    # 4
    for v in eord.keys():
        random.shuffle(eord[v])

    # 5
    # a remove arc of eord in T
    for i in eord.keys():
        if i != lst:
            eord[i].remove(T[i][0])
    # b
    vertex = fst
    seq_perm = vertex
    step = 0
    while step <= len_seq - k:  # change into raise error
        vertex = random_seq(eord, vertex, T)
        seq_perm += vertex[-1]
        step += 1

    assert same_klets(seq_perm, seq, k)
    print(seq_perm)


def main():
    S1 = "GTTGCGAAGCCTACACTGATATATGAATCCAAGCTAGAGCAGGGCTCTTAAAATTCGGAGTTGTAGATGCTCAATACTCCAATCGGTTTTTTCGTGCGTT"
    S1 = "ATCAGCATGCCGTATC"
    print(S1)
    victor_algorithm(S1, 2)
    victor_algorithm(S1, 3)
    victor_algorithm(S1, 4)
    # test if uniform
    # check T
    # test linearity of the algo


if __name__ == "__main__":
    main()
