#!/usr/bin/env python3
"""Generate a synthetic VCF file with 100,000 variants on chr22 for benchmarking."""

import random

random.seed(42)

NUM_VARIANTS = 100_000
CHROM = "chr22"
POS_START = 10_000_000
POS_END = 51_000_000
OUTPUT = "/Users/davidcho/Projects/VEP_ensemble_C/test_data/e2e/benchmark_100k.vcf"

BASES = ["A", "C", "G", "T"]

def random_alt_snv(ref):
    """Return a random base different from ref."""
    alts = [b for b in BASES if b != ref]
    return random.choice(alts)

def random_insertion():
    """Return a random insertion sequence (1-10 bases)."""
    length = random.randint(1, 10)
    return "".join(random.choice(BASES) for _ in range(length))

def random_deletion_length():
    """Return a random deletion length (1-10 bases)."""
    return random.randint(1, 10)

def main():
    # Generate sorted unique positions
    positions = sorted(random.sample(range(POS_START, POS_END + 1), NUM_VARIANTS))

    with open(OUTPUT, "w") as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr22>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for pos in positions:
            variant_type = random.random()

            if variant_type < 0.70:
                # SNV (70%)
                ref = random.choice(BASES)
                alt = random_alt_snv(ref)
            elif variant_type < 0.85:
                # Insertion (15%): anchor base + inserted sequence
                ref = random.choice(BASES)
                alt = ref + random_insertion()
            else:
                # Deletion (15%): multi-base ref, single-base alt (anchor)
                del_len = random_deletion_length()
                ref = "".join(random.choice(BASES) for _ in range(1 + del_len))
                alt = ref[0]

            f.write(f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"Wrote {NUM_VARIANTS} variants to {OUTPUT}")

if __name__ == "__main__":
    main()
