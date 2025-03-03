#!/usr/bin/env python3

import sys
import pickle

def compare_suffix_naive(genome, pattern, suffix_start):
    """
    Compare 'pattern' to the substring genome[suffix_start:] using naive character-by-character
    comparison. Return:
      -1 if pattern < suffix, 0 if pattern == suffix (prefix-match), +1 if pattern > suffix
      char_used is how many character comparisons were made.
    """
    i = 0
    max_len = min(len(pattern), len(genome) - suffix_start)
    while i < max_len:
        if pattern[i] < genome[suffix_start + i]:
            return -1, i+1  # we used (i+1) comparisons
        elif pattern[i] > genome[suffix_start + i]:
            return +1, i+1
        i += 1

    # If we exhaust the entire pattern or suffix portion:
    if len(pattern) == max_len:
        # Pattern length == max_len
        if max_len == (len(genome) - suffix_start):
            return 0, i  # pattern and suffix are exactly equal length
        else:
            return -1, i  # pattern is prefix of suffix, so pattern < suffix
    else:
        # Suffix portion ended first
        return +1, i


def compare_suffix_simpaccel(genome, pattern, suffix_start, offset):
    """
    A 'simple accelerant' version:
      - We start comparing from `offset` (the minimum of LCP with LB/UB).
      - Return the same tuple as compare_suffix_naive, except we skip the first `offset` chars.
    """
    i = offset
    max_len = min(len(pattern), len(genome) - suffix_start)
    comparisons_used = 0  # track how many char comparisons we do *this call*

    while i < max_len:
        comparisons_used += 1
        if pattern[i] < genome[suffix_start + i]:
            return -1, comparisons_used
        elif pattern[i] > genome[suffix_start + i]:
            return +1, comparisons_used
        i += 1

    # If we exhaust the entire pattern or suffix portion:
    if len(pattern) == max_len:
        if max_len == (len(genome) - suffix_start):
            return 0, comparisons_used  # exact match
        else:
            return -1, comparisons_used  # pattern is prefix of suffix
    else:
        return +1, comparisons_used


def binary_search_bound(genome, sa, pattern, mode="naive", find_lower_bound=True):
    """
    Performs a standard binary search to find either the lower bound or upper bound
    of 'pattern' in suffix array 'sa'. Returns (bound_index, total_char_comps).

    If find_lower_bound=True, we find the first suffix >= pattern.
    If find_lower_bound=False, we find the first suffix > pattern.

    'mode' can be 'naive' or 'simpaccel'.
    """
    left, right = 0, len(sa)
    char_comps = 0
    # We maintain "lcpLB" and "lcpUB" for simpaccel
    lcp_left = 0
    lcp_right = 0

    while left < right:
        mid = (left + right) // 2
        suffix_start = sa[mid]

        # Choose which compare function to call
        if mode == "naive":
            result, used = compare_suffix_naive(genome, pattern, suffix_start)
        else:
            # We start comparing at offset = min(lcp_left, lcp_right)
            offset = min(lcp_left, lcp_right)
            result, used = compare_suffix_simpaccel(genome, pattern, suffix_start, offset)

        char_comps += used

        if result < 0 or (result == 0 and find_lower_bound):
            # suffix >= pattern => go left if we want the lower bound
            right = mid
            # Update lcp_right if using simpaccel
            if mode == "simpaccel":
                # lcp_right should be offset + number of matched chars
                # but if we broke early due to mismatch, matched chars = used - 1
                # We'll do a naive approach: if result < 0 => we matched 'used-1' chars
                # if result == 0 => matched 'used' chars
                # For simplicity, let's estimate:
                # We'll treat "matched" as (offset + used) if result == 0,
                # else offset + (used - 1).  This is a rough approach.
                if result == 0:
                    matched = offset + used
                else:
                    matched = offset + (used - 1)
                if matched < 0:
                    matched = 0
                lcp_right = max(0, matched)
        else:
            # suffix < pattern => go right
            left = mid + 1
            if mode == "simpaccel":
                # similarly update lcp_left
                if result == 0:
                    matched = offset + used
                else:
                    matched = offset + (used - 1)
                if matched < 0:
                    matched = 0
                lcp_left = max(0, matched)

    return left, char_comps


def query(genome, sa, pattern, mode="naive"):
    """
    Returns (lb, ub, lb_cmp, ub_cmp) for the pattern:
      - lb, ub: the range [lb, ub) of suffix array indices where 'pattern' occurs.
      - lb_cmp, ub_cmp: how many character comparisons were done for each binary search.
    We do two searches with pattern# and pattern}.
    """
    # form P# and P}
    lower_pat = pattern + "#"
    upper_pat = pattern + "}"

    lb, lb_cmp = binary_search_bound(genome, sa, lower_pat, mode, find_lower_bound=True)
    ub, ub_cmp = binary_search_bound(genome, sa, upper_pat, mode, find_lower_bound=True)
    return lb, ub, lb_cmp, ub_cmp


def parse_fasta(fasta_path):
    """
    Yields (header, sequence) for each record in the FASTA.
    Handles multi-line sequences until next '>' line.
    """
    header = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # output previous record
                if header is not None:
                    yield (header, "".join(seq_chunks))
                header = line[1:]  # remove '>'
                seq_chunks = []
            else:
                # accumulate sequence lines
                seq_chunks.append(line)
        # end of file => yield last record if any
        if header is not None:
            yield (header, "".join(seq_chunks))


def query_sa(index_path, queries_fasta, query_mode, output_path):
    """
    Main function:
      1) Load genome, sa from index_path
      2) Parse each record from queries_fasta
      3) For each record, call query(...)
      4) Write to output_path in the required format
    """
    with open(index_path, 'rb') as pf:
        genome, sa = pickle.load(pf)

    with open(output_path, 'w') as ofile:
        for (qname, qseq) in parse_fasta(queries_fasta):
            lb, ub, lb_cmp, ub_cmp = query(genome, sa, qseq, mode=query_mode)
            k = ub - lb  # number of occurrences
            ofile.write(f"{qname}\t{lb_cmp}\t{ub_cmp}\t{k}")
            if k > 0:
                hits = [str(sa[i]) for i in range(lb, ub)]
                ofile.write("\t" + "\t".join(hits))
            ofile.write("\n")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python querysa.py <index_path> <queries_fasta> <query_mode> <output_path>")
        sys.exit(1)

    index_path = sys.argv[1]
    queries_fasta = sys.argv[2]
    query_mode = sys.argv[3]
    if query_mode not in ("naive", "simpaccel"):
        print("Error: query_mode must be 'naive' or 'simpaccel'.")
        sys.exit(1)
    output_path = sys.argv[4]

    query_sa(index_path, queries_fasta, query_mode, output_path)
