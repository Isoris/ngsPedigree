#!/usr/bin/env python3
"""
Generate a synthetic ngsRelate .res + samples.txt fixture with hand-checked
truth, exercising every topology branch the classifier needs to handle.

Design (12 samples, 5 hubs):

  Hub H001: two_parents_with_sibship
    P1 (CGA001) and P2 (CGA002) — unrelated to each other (theta ~ 0)
    Offspring: CGA003, CGA004, CGA005 (full sibs)
    P1-{O1,O2,O3}: PO   P2-{O1,O2,O3}: PO   O-O: FS

  Hub H002: parent_with_sibship (one parent only)
    P (CGA006)
    Offspring: CGA007, CGA008, CGA009
    All 3 offspring full sibs to each other

  Hub H003: sibship_only (parents not sampled)
    CGA010, CGA011 — full sibs to each other

  Hub H004: po_dyad_only (single PO edge, no triangulation)
    CGA012 — paired with CGA001 as PO. Wait, that would put CGA012 in H001.
    Use a separate dyad instead:
    Actually with only 12 samples we'll skip po_dyad_only here and add it
    in a second smaller fixture. Let's keep this fixture clean.

  Plus: duplicate pair would need 13 samples. Same issue — keep clean.

  Actually, let me redesign: use 13 samples with a duplicate, and document
  po_dyad_only as a tested branch via a second tiny fixture.

Final design (13 samples, 4 hubs + 1 isolated):

  Hub H001 (two_parents_with_sibship): CGA001, CGA002 + CGA003, CGA004, CGA005
  Hub H002 (parent_with_sibship):      CGA006 + CGA007, CGA008, CGA009
  Hub H003 (sibship_only):             CGA010, CGA011
  Hub H004 (duplicate_pair):           CGA012, CGA013
  Isolated: none — all 13 in some hub.

  Plus a separate po_dyad fixture with 2 samples: CGA014, CGA015.
"""

import sys
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent

# Sample IDs in ngsRelate index order.
SAMPLES = [f"CGA{i:03d}" for i in range(1, 14)]  # CGA001..CGA013

# (sample_a, sample_b, edge_truth) where edge_truth is what we expect
# the classifier to produce.
TRUTH_EDGES = {
    # H001: CGA001, CGA002 unrelated parents; CGA003-005 offspring of both
    ("CGA001", "CGA002"): "unrelated",
    ("CGA001", "CGA003"): "parent_offspring",
    ("CGA001", "CGA004"): "parent_offspring",
    ("CGA001", "CGA005"): "parent_offspring",
    ("CGA002", "CGA003"): "parent_offspring",
    ("CGA002", "CGA004"): "parent_offspring",
    ("CGA002", "CGA005"): "parent_offspring",
    ("CGA003", "CGA004"): "full_sibling",
    ("CGA003", "CGA005"): "full_sibling",
    ("CGA004", "CGA005"): "full_sibling",
    # H002: CGA006 parent; CGA007-009 offspring
    ("CGA006", "CGA007"): "parent_offspring",
    ("CGA006", "CGA008"): "parent_offspring",
    ("CGA006", "CGA009"): "parent_offspring",
    ("CGA007", "CGA008"): "full_sibling",
    ("CGA007", "CGA009"): "full_sibling",
    ("CGA008", "CGA009"): "full_sibling",
    # H003: pure sibship
    ("CGA010", "CGA011"): "full_sibling",
    # H004: duplicate
    ("CGA012", "CGA013"): "duplicate_or_clone",
}


def synth_row(sample_a, sample_b, edge_truth, samples):
    """Generate a synthetic ngsRelate row with values consistent with edge_truth."""
    a_idx = samples.index(sample_a)
    b_idx = samples.index(sample_b)

    # Defaults for unrelated:
    theta = 0.005
    ibs0 = 0.18
    king = 0.0
    r0 = 0.5
    r1 = 0.5
    j9 = 0.05
    n_sites = 50000

    if edge_truth == "parent_offspring":
        theta = 0.247
        ibs0 = 0.0008
        king = 0.245
        r0 = 0.02
        r1 = 0.85
        j9 = 0.05
    elif edge_truth == "full_sibling":
        theta = 0.255
        ibs0 = 0.045
        king = 0.250
        r0 = 0.18
        r1 = 0.62
        j9 = 0.22
    elif edge_truth == "duplicate_or_clone":
        theta = 0.495
        ibs0 = 0.0001
        king = 0.495
        r0 = 0.01
        r1 = 0.01
        j9 = 0.97
    elif edge_truth == "second_degree":
        theta = 0.12
        ibs0 = 0.10
        king = 0.12
        r0 = 0.32
        r1 = 0.45
        j9 = 0.10
    elif edge_truth == "third_degree":
        theta = 0.06
        ibs0 = 0.14
        king = 0.06
        r0 = 0.42
        r1 = 0.32
        j9 = 0.05
    # else: unrelated keeps defaults.

    return dict(a=a_idx, b=b_idx, nSites=n_sites,
                J7=0.05, J8=0.05, J9=j9,
                rab=2*theta, Fa=0.01, Fb=0.01,
                theta=theta, inbreed_a=0.01, inbreed_b=0.01,
                **{"2of3_IDB": 0.5}, FDiff=0.0,
                loglh=-1000.0, nIter=20, coverage=9.0,
                IBS0=ibs0, IBS1=0.5, IBS2=0.5-ibs0,
                R0=r0, R1=r1, KING=king)


def write_samples(path):
    with open(path, "w") as fh:
        for s in SAMPLES:
            fh.write(s + "\n")


def write_res(path):
    cols = ["a", "b", "nSites", "J7", "J8", "J9", "rab", "Fa", "Fb",
            "theta", "inbreed_a", "inbreed_b", "2of3_IDB", "FDiff",
            "loglh", "nIter", "coverage", "IBS0", "IBS1", "IBS2",
            "R0", "R1", "KING"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, sa in enumerate(SAMPLES):
            for j, sb in enumerate(SAMPLES):
                if j <= i:
                    continue
                edge_truth = TRUTH_EDGES.get((sa, sb)) or TRUTH_EDGES.get((sb, sa)) or "unrelated"
                row = synth_row(sa, sb, edge_truth, SAMPLES)
                vals = [str(row[c]) for c in cols]
                fh.write("\t".join(vals) + "\n")


def write_truth(path):
    """Write the expected hub roster for verification."""
    truth = {
        # H001
        "CGA001": ("two_parents_with_sibship", "parent_a"),
        "CGA002": ("two_parents_with_sibship", "parent_b"),
        "CGA003": ("two_parents_with_sibship", "possible_offspring"),
        "CGA004": ("two_parents_with_sibship", "possible_offspring"),
        "CGA005": ("two_parents_with_sibship", "possible_offspring"),
        # H002
        "CGA006": ("parent_with_sibship", "forced_parent"),
        "CGA007": ("parent_with_sibship", "possible_offspring"),
        "CGA008": ("parent_with_sibship", "possible_offspring"),
        "CGA009": ("parent_with_sibship", "possible_offspring"),
        # H003
        "CGA010": ("sibship_only", "possible_full_sib"),
        "CGA011": ("sibship_only", "possible_full_sib"),
        # H004
        "CGA012": ("duplicate_pair", "duplicate"),
        "CGA013": ("duplicate_pair", "duplicate"),
    }
    with open(path, "w") as fh:
        fh.write("sample_id\texpected_hub_type\texpected_role\n")
        for sid, (ht, role) in truth.items():
            fh.write(f"{sid}\t{ht}\t{role}\n")


def write_sex(path):
    """Write a sex TSV that should promote some roles."""
    sex = {
        "CGA001": "male",     # → father in H001
        "CGA002": "female",   # → mother in H001
        "CGA006": "female",   # → mother in H002
        # CGA003..009: leave unspecified (offspring; sex doesn't change role)
        # CGA010..013: leave unspecified
    }
    with open(path, "w") as fh:
        fh.write("sample_id\tsex\n")
        for sid, sx in sex.items():
            fh.write(f"{sid}\t{sx}\n")


def main():
    write_samples(OUT_DIR / "samples.txt")
    write_res(OUT_DIR / "synth.res")
    write_truth(OUT_DIR / "expected_truth.tsv")
    write_sex(OUT_DIR / "sex.tsv")
    print(f"Wrote synthetic fixture to {OUT_DIR}")
    print(f"  samples.txt        ({len(SAMPLES)} samples)")
    print(f"  synth.res          ({len(SAMPLES) * (len(SAMPLES) - 1) // 2} pairs)")
    print(f"  expected_truth.tsv (per-sample expected hub_type and role)")
    print(f"  sex.tsv            (3 entries; should promote 3 parent labels)")


if __name__ == "__main__":
    main()
