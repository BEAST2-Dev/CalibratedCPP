#!/usr/bin/env python3
"""
Compute the Exponential mean for each of the paper's offset-exponential calibrations, given
the offset (= hard minimum bound) and soft maximum bound from de Vries & Beck Table 1.

The paper states its convention explicitly: "the shape of the exponential distribution
[is] specified such that there is a 5% probability of the divergence being older than the
maximum bound." For an offset-exponential age = offset + Exp(mean), the survival function
of Exp(mean) is P(X > x) = exp(-x/mean). Setting P(age > upper) = 0.05 and solving for mean:

    exp(-(upper - offset) / mean) = 0.05
    mean = (upper - offset) / -ln(0.05)
         = (upper - offset) / 2.995732...

Usage:
    python compute_offset_exponential_mean.py
"""
import math

TAIL_PROBABILITY = 0.05  # P(age > upper_bound), per the paper's stated convention
DENOM = -math.log(TAIL_PROBABILITY)  # ~2.995732

# (node number, clade, offset [= hard lower bound], upper [= soft maximum bound]), Ma.
# The 6 offset-exponential calibrations present in this project's 29-taxon dataset (Glires,
# the 7th node the paper lists this shape for, isn't sampled here).
NODES = [
    (1,  "Euarchontoglires", 65.79,  125.816),
    (3,  "Euarchonta",       65.79,  125.816),
    (5,  "Primates",         55.935, 66.095),
    (13, "Cercopithecidae",  12.47,  25.235),
    (18, "Hominoidea",       13.4,   25.235),
    (19, "Hominidae",        12.3,   25.235),
]


def compute_mean(offset, upper, tail_probability=TAIL_PROBABILITY):
    return (upper - offset) / -math.log(tail_probability)


def main():
    print(f"Tail probability P(age > upper) = {TAIL_PROBABILITY}  ->  denom = -ln({TAIL_PROBABILITY}) = {DENOM:.6f}\n")
    print(f"{'Node':>5} {'Clade':<20} {'offset':>10} {'upper':>10} {'mean':>12}")
    print("-" * 62)
    for node, clade, offset, upper in NODES:
        mean = compute_mean(offset, upper)
        print(f"{node:>5} {clade:<20} {offset:>10} {upper:>10} {mean:>12.6f}")


if __name__ == "__main__":
    main()
