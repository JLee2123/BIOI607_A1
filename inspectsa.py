print("DEBUG: In inspectsa.py, before anything else!")

def inspect_sa(index_path, sample_rate, output_path):
    import pickle

    print(f"Attempting to load index from {index_path}...")
    (genome, sa) = pickle.load(open(index_path, 'rb'))
    print(f"Loaded genome of length {len(genome)} and suffix array of length {len(sa)}")

    # Edge case: if there's only one suffix or empty SA, LCP statistics are trivially zero.
    if len(sa) < 2:
        mean_lcp = 0.0
        median_lcp = 0.0
        max_lcp = 0
        indices_to_sample = [0] if len(sa) == 1 else []
        spot_check_values = [str(sa[i]) for i in indices_to_sample]
    else:
        def compute_lcp(g, s1, s2):
            i = 0
            while s1 + i < len(g) and s2 + i < len(g) and g[s1 + i] == g[s2 + i]:
                i += 1
            return i

        # Compute LCP for adjacent suffixes
        lcp1_vals = []
        for i in range(len(sa) - 1):
            lcp_val = compute_lcp(genome, sa[i], sa[i+1])
            lcp1_vals.append(lcp_val)

        # mean LCP
        mean_lcp = sum(lcp1_vals) / len(lcp1_vals)

        # median LCP
        lcp1_vals_sorted = sorted(lcp1_vals)
        mid = len(lcp1_vals_sorted) // 2
        if len(lcp1_vals_sorted) % 2 == 1:
            median_lcp = float(lcp1_vals_sorted[mid])
        else:
            median_lcp = (lcp1_vals_sorted[mid - 1] + lcp1_vals_sorted[mid]) / 2.0

        # max LCP
        max_lcp = max(lcp1_vals)

        # spot check sampling
        limit = len(sa) - 1
        max_k = len(sa) // sample_rate
        indices_to_sample = [sample_rate * i for i in range(max_k + 1)]
        indices_to_sample = [idx for idx in indices_to_sample if idx <= limit]
        spot_check_values = [str(sa[idx]) for idx in indices_to_sample]

    print(f"Writing statistics to {output_path}...")
    with open(output_path, 'w') as ofile:
        ofile.write(f"{mean_lcp}\n")      # 1) mean LCP
        ofile.write(f"{median_lcp}\n")    # 2) median LCP
        ofile.write(f"{max_lcp}\n")       # 3) max LCP
        ofile.write("\t".join(spot_check_values) + "\n")  # 4) spot checks

    print(f"Done! Check {output_path} for output.")


if __name__ == "__main__":
    import sys
    print("DEBUG: In the __main__ block!")

    # Ensure the user provided three arguments: index_path, sample_rate, and output_path
    if len(sys.argv) < 4:
        print("Usage: python my_inspectsa.py <index_path> <sample_rate> <output_path>")
        sys.exit(1)

    index_path = sys.argv[1]
    sample_rate = int(sys.argv[2])
    output_path = sys.argv[3]

    # Call the function
    inspect_sa(index_path, sample_rate, output_path)
