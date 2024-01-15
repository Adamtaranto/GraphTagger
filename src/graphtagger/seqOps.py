def revComp(seq):
    """Rev comp DNA string."""
    revcompl = lambda x: "".join(
        [{"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}[B] for B in x][::-1]
    )
    return revcompl(seq)
