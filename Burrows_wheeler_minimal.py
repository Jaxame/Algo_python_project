'''
 # @ Author: Claire Lemaître & Pierre Peterlongo
 # @ Creation: 2021-11-21 
 # @ Last update: 2021-04-01
 # @ Python Version : 3
 '''

from typing import Dict, List, Tuple


def get_bwt(s: str, sa: List[int]) -> str:
    """
    Returns the Burrows-Wheeler Transform (BWT) of a sequence.

    Args:
        s (str): [description]
        sa (List[int]): [description]

    Returns:
        bwt(str): [description]
    """
    bwt = ''
    for pos in sa:
        if pos == 0:
            bwt += '$'
        else:
            bwt += s[pos - 1]
    return bwt


def get_n(bwt: str) -> Dict[str, int]:
    """
    Returns the number of occurrences of each character in BWT.

    Args:
        bwt (str): (any) sequence from which we compute n

    Returns:
        Dict[str, int]: the number of occurrences of each 
        character in bwt
    """
    # Dictionnary comprehension (directly in constructor)
    n = {letter: 0 for letter in 'ATGC$'}
    for letter in bwt:
        n[letter] += 1
    return n


def get_ranks(bwt: str) -> List[int]:
    """
    Returns the rank of each character in the given BWT.

    Args:
        bwt (str): (any) sequence from which we compute in

    Returns:
        List[int]: rank of each character in bwt
    """
    r = []
    n = {letter: 0 for letter in 'ATGC$'}
    for letter in bwt:
        n[letter] += 1
        r.append(n[letter])
    return r


def left_first(alpha: str, k: int, n: Dict[str, int]) -> int:
    """
    Returns the line of the k^th occurrence of character alpha in the F.

    Args:
        alpha (str): concerned character
        k (int): k^th occurrence of alpha
        n (dict): line of the first occurrence of alpha

    Raises:
        KeyError: if alpha not in alphabet

    Returns:
        int: line corresponding to the k^th occurrence
             of character alpha in F
    """
    if alpha == '$':
        return 0
    if alpha == 'A':
        return n['$'] + k - 1
    if alpha == 'C':
        return n['$'] + n['A'] + k - 1
    if alpha == 'G':
        return n['$'] + n['A'] + n['C'] + k - 1
    if alpha == 'T':
        return n['$'] + n['A'] + n['C'] + n['G'] + k - 1
    raise KeyError


def bwt_2_seq(bwt: str, n: Dict[str, int], r: List[int]) -> str:
    """
    Reverses the BWT to the original sequence.

    Args:
        bwt (str): BWT
        n (dict): line of the first occurrence of alpha
        r (list of int): rank of each character in bwt

    Returns:
        str: original sequence
    """
    sequence = ''
    line = 0
    alpha = bwt[line]
    while alpha != '$':
        sequence = alpha + sequence
        line = left_first(alpha, r[line], n)
        alpha = bwt[line]
    return sequence


def first_line_in_bwt_containing_alpha(bwt: str, start: int, max_line: int,
                                       alpha: str) -> int:
    """
    Returns first BWT line index containing alpha.

    Args:
        bwt (str): BWT
        start (int): start index
        max_line (int): max index
        alpha (str): character to look for

    Returns:
        int: first BWT line index containing alpha, else max_line
    """
    line = start
    while line < max_line + 1 and bwt[line] != alpha:
        line += 1
    return line


def last_line_in_bwt_containing_alpha(bwt: str, stop: int, min_line: int,
                                      alpha: str) -> str:
    """
    Returns the last BWT line index containing alpha.

    Args:
        bwt (str): BWT
        stop (int): stop index
        min_line (int): min index
        alpha (str): character to look for

    Returns:
        int: last BWT line index containing alpha, else min_line
    """
    line = stop
    while line > min_line and bwt[line] != alpha:
        line -= 1
    return line


def contains(p: str, bwt: str, n: Dict[str, int],
             r: List[int]) -> Tuple[bool, int, int]:
    """
    Returns True if pattern p exists in sequence s.

    Args:
        p (str): pattern to find
        bwt (str): BWT
        n (dict): dictionary of occurrences
        r (list): list of ranks in the BWT

    Returns:
        tuple: [description]
    """
    start = 0
    stop = len(bwt) - 1

    for i in range(len(p) - 1, -1, -1):
        # Finding the interval's beginning
        new_start = first_line_in_bwt_containing_alpha(bwt, start, stop, p[i])
        if new_start > stop:  # game over, p not in s
            return False, -1, -1

        # Finding the interval's end
        new_stop = last_line_in_bwt_containing_alpha(bwt, stop, start, p[i])
        start = left_first(p[i], r[new_start], n)
        stop = left_first(p[i], r[new_stop], n)
    return True, start, stop