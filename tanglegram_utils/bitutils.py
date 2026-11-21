"""Helpers for leaf sets represented as integers (bitmasks)."""

# popcount (int.bit_count exists since Python 3.10)
pc = int.bit_count if hasattr(int, "bit_count") else (lambda x: bin(x).count("1"))


def bits(x):
    while x:
        lb = x & -x
        yield lb.bit_length() - 1
        x ^= lb
