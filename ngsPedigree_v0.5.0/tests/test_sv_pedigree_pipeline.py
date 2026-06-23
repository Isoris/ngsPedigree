"""
test_sv_pedigree_pipeline — end-to-end test for the broke-grad-student
pipeline: Delly + Manta DEL VCFs → merged catalogue → KING-robust
relatedness → Mendelian exclusion-based parent ID → chromosome
inheritance map.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.catalogue_merge import (  # noqa: E402
    cohort_samples,
    merge_two_callers,
    records_match,
    to_genotype_matrix,
)
from hpp.del_inheritance import (  # noqa: E402
    DelMarkerLocus,
    build_inheritance_map_for_dyad,
    build_inheritance_map_for_triad,
    exclude_as_parent,
    transmitted_from_parent,
    transmitted_from_parent_triad,
)
from hpp.del_relatedness import all_pairs, pair_coeffs  # noqa: E402
from hpp.relatedness_sim import classify_edge_stdlib  # noqa: E402
from hpp.vcf_sv import normalise_gt, read_del_calls  # noqa: E402

FIX = THIS_DIR / "fixtures" / "synthetic_svvcfs"


# ----------------------------------------------------------------------
# VCF SV reader.
# ----------------------------------------------------------------------


class TestVcfReader(unittest.TestCase):
    def test_delly_pass_only(self):
        recs = read_del_calls(FIX / "delly.vcf")
        # 16 records in fixture, 1 LowQual → 15 PASS DEL records
        self.assertEqual(len(recs), 15)
        self.assertTrue(all(r.sv_type == "DEL" for r in recs))
        self.assertTrue(all(r.filter_pass for r in recs))

    def test_manta_records_counted(self):
        recs = read_del_calls(FIX / "manta.vcf")
        self.assertEqual(len(recs), 9)

    def test_caller_detected(self):
        d = read_del_calls(FIX / "delly.vcf")
        m = read_del_calls(FIX / "manta.vcf")
        self.assertTrue(all(r.caller == "delly" for r in d))
        self.assertTrue(all(r.caller == "manta" for r in m))

    def test_genotypes_parsed(self):
        recs = read_del_calls(FIX / "delly.vcf")
        first = recs[0]
        self.assertEqual(first.chrom, "Chr1")
        self.assertEqual(first.pos, 1000)
        self.assertEqual(first.end, 1500)
        self.assertEqual(first.genotypes["P_F"], "0/1")
        self.assertEqual(first.genotypes["P_M"], "0/0")

    def test_gt_normalises_phased(self):
        self.assertEqual(normalise_gt("0|1"), "0/1")
        self.assertEqual(normalise_gt("1|1"), "1/1")
        self.assertEqual(normalise_gt("."), "./.")
        self.assertEqual(normalise_gt(None), "./.")


# ----------------------------------------------------------------------
# Catalogue merger.
# ----------------------------------------------------------------------


class TestCatalogueMerge(unittest.TestCase):
    def setUp(self):
        self.delly = read_del_calls(FIX / "delly.vcf")
        self.manta = read_del_calls(FIX / "manta.vcf")

    def test_matching_records_found(self):
        # Delly DEL @ Chr1:1000-1500 should match Manta @ Chr1:1050-1530
        d = next(r for r in self.delly if r.chrom == "Chr1" and r.pos == 1000)
        m = next(r for r in self.manta if r.chrom == "Chr1" and r.pos == 1050)
        self.assertTrue(records_match(d, m, bp_tolerance=500,
                                       reciprocal_overlap=0.5))

    def test_merged_catalogue_size(self):
        merged = merge_two_callers(self.delly, self.manta)
        # Delly: 15; Manta: 9. Overlapping pairs: at least the explicit
        # matches at (1000/1050), (10000/10010), (20000/20020),
        # (50000/50010), (5000/5020), (10000/10005), (50000-Chr2 missing),
        # (10000/10005). Exact merged count depends on the matcher.
        self.assertGreaterEqual(len(merged), 15)
        self.assertLessEqual(len(merged), 24)

    def test_merged_catalogue_has_both_caller_markers(self):
        merged = merge_two_callers(self.delly, self.manta)
        # Should find at least one marker with both callers.
        n_both = sum(1 for m in merged if m.n_callers == 2)
        self.assertGreater(n_both, 0)

    def test_cohort_samples_recovered(self):
        merged = merge_two_callers(self.delly, self.manta)
        samples = cohort_samples(merged)
        self.assertEqual(set(samples), {"P_F", "P_M", "C", "UNREL"})


# ----------------------------------------------------------------------
# KING-robust relatedness.
# ----------------------------------------------------------------------


class TestKingRobust(unittest.TestCase):
    def setUp(self):
        delly = read_del_calls(FIX / "delly.vcf")
        manta = read_del_calls(FIX / "manta.vcf")
        merged = merge_two_callers(delly, manta)
        self.gmatrix = to_genotype_matrix(merged)
        self.samples = cohort_samples(merged)

    def test_unrelated_pair_has_low_theta(self):
        # UNREL was constructed to be Mendelian-independent from everyone.
        p = pair_coeffs("P_F", "UNREL", self.gmatrix)
        self.assertIsNotNone(p.theta)
        # KING-robust returns NEGATIVE values for unrelated pairs.
        self.assertLess(p.theta, 0.1)

    def test_parent_offspring_pair_above_first_degree(self):
        # P_F → C and P_M → C should both be classified as first-degree.
        p = pair_coeffs("P_F", "C", self.gmatrix)
        self.assertIsNotNone(p.theta)
        # KING-robust kinship for true PO is ≈ 0.25 in expectation;
        # with our small marker set it lands above 0.15.
        self.assertGreater(p.theta, 0.10)

    def test_classifier_calls_PO(self):
        p = pair_coeffs("P_F", "C", self.gmatrix)
        v = classify_edge_stdlib({
            "theta": p.theta, "IBS0": p.IBS0,
            "KING": p.theta, "nSites": p.n_informative,
        })
        # First-degree by theta + low IBS0 → parent_offspring
        self.assertIn(v["edge_class"], {"parent_offspring", "full_sibling"})

    def test_unrelated_classified_as_unrelated(self):
        p = pair_coeffs("P_F", "UNREL", self.gmatrix)
        v = classify_edge_stdlib({
            "theta": p.theta, "IBS0": p.IBS0,
            "KING": p.theta, "nSites": p.n_informative,
        })
        self.assertIn(v["edge_class"], {"unrelated", "third_degree"})


# ----------------------------------------------------------------------
# Mendelian exclusion-based parent identification.
# ----------------------------------------------------------------------


class TestExclusion(unittest.TestCase):
    def setUp(self):
        merged = merge_two_callers(
            read_del_calls(FIX / "delly.vcf"),
            read_del_calls(FIX / "manta.vcf"),
        )
        self.gmatrix = to_genotype_matrix(merged)

    def test_true_parents_pass_exclusion(self):
        # P_F is the true parent of C — no opposite-hom contradictions.
        e = exclude_as_parent("P_F", "C", self.gmatrix)
        self.assertEqual(e.n_excluding_markers, 0)
        self.assertTrue(e.can_be_parent)
        e2 = exclude_as_parent("P_M", "C", self.gmatrix)
        self.assertEqual(e2.n_excluding_markers, 0)
        self.assertTrue(e2.can_be_parent)

    def test_unrelated_fails_exclusion(self):
        # UNREL was deliberately set hom_alt at a marker where C is hom_ref.
        # That makes UNREL incompatible as C's parent.
        e = exclude_as_parent("UNREL", "C", self.gmatrix)
        self.assertGreater(e.n_excluding_markers, 0)
        self.assertFalse(e.can_be_parent)


# ----------------------------------------------------------------------
# Transmitted-allele inference.
# ----------------------------------------------------------------------


class TestTransmission(unittest.TestCase):
    def test_hom_parent_transmits_its_allele(self):
        self.assertEqual(transmitted_from_parent("0/0", "0/0"), "REF")
        self.assertEqual(transmitted_from_parent("0/0", "0/1"), "REF")
        self.assertEqual(transmitted_from_parent("1/1", "1/1"), "DEL")
        self.assertEqual(transmitted_from_parent("1/1", "0/1"), "DEL")

    def test_het_parent_to_homozygous_offspring_resolvable(self):
        self.assertEqual(transmitted_from_parent("0/1", "0/0"), "REF")
        self.assertEqual(transmitted_from_parent("0/1", "1/1"), "DEL")

    def test_het_parent_het_offspring_ambiguous(self):
        self.assertEqual(transmitted_from_parent("0/1", "0/1"), "ambiguous")

    def test_contradictions(self):
        self.assertEqual(transmitted_from_parent("0/0", "1/1"), "contradiction")
        self.assertEqual(transmitted_from_parent("1/1", "0/0"), "contradiction")

    def test_triad_disambiguates_het_het(self):
        # parent HET, offspring HET, co-parent HOM_REF → parent must have given DEL
        self.assertEqual(
            transmitted_from_parent_triad("0/1", "0/0", "0/1"), "DEL")
        # parent HET, offspring HET, co-parent HOM_INV → parent gave REF
        self.assertEqual(
            transmitted_from_parent_triad("0/1", "1/1", "0/1"), "REF")
        # co-parent also HET → still ambiguous
        self.assertEqual(
            transmitted_from_parent_triad("0/1", "0/1", "0/1"), "ambiguous")


# ----------------------------------------------------------------------
# Chromosome inheritance map.
# ----------------------------------------------------------------------


class TestInheritanceMap(unittest.TestCase):
    def setUp(self):
        merged = merge_two_callers(
            read_del_calls(FIX / "delly.vcf"),
            read_del_calls(FIX / "manta.vcf"),
        )
        self.gmatrix = to_genotype_matrix(merged)
        self.loci = [
            DelMarkerLocus(m.marker_id, m.chrom, m.start, m.end)
            for m in merged
        ]

    def test_segments_emitted_for_dyad(self):
        segs = build_inheritance_map_for_dyad(
            dyad_id="dyad_PF_C",
            parent_sample_id="P_F",
            offspring_sample_id="C",
            loci=self.loci,
            genotype_matrix=self.gmatrix,
        )
        # At least one segment per chromosome in the catalogue.
        self.assertGreater(len(segs), 0)
        chroms = {s.chrom for s in segs}
        self.assertEqual(chroms, {"Chr1", "Chr2"})

    def test_segment_fields_well_formed(self):
        segs = build_inheritance_map_for_dyad(
            dyad_id="dyad_PF_C",
            parent_sample_id="P_F",
            offspring_sample_id="C",
            loci=self.loci,
            genotype_matrix=self.gmatrix,
        )
        for s in segs:
            self.assertIn(s.transmitted_allele,
                          {"REF", "DEL", "ambiguous"})
            self.assertIn(s.confidence, {"Gold", "Silver", "Bronze"})
            self.assertIn(s.parental_hap_inherited,
                          {"1", "2", "ambiguous"})
            self.assertLessEqual(s.seg_start, s.seg_end)

    def test_triad_map_two_parental_traces(self):
        m = build_inheritance_map_for_triad(
            triad_id="triad_T",
            paternal_sample_id="P_F",
            maternal_sample_id="P_M",
            offspring_sample_id="C",
            loci=self.loci,
            genotype_matrix=self.gmatrix,
        )
        self.assertIn("paternal", m)
        self.assertIn("maternal", m)
        self.assertGreater(len(m["paternal"]), 0)
        self.assertGreater(len(m["maternal"]), 0)


# ----------------------------------------------------------------------
# End-to-end CLI smoke.
# ----------------------------------------------------------------------


class TestCliEndToEnd(unittest.TestCase):
    def test_full_pipeline_runs(self):
        script = THIS_DIR.parent / "scripts" / "08_pedigree_from_sv_vcfs.py"
        with tempfile.TemporaryDirectory() as td:
            outdir = Path(td)
            res = subprocess.run(
                [sys.executable, str(script),
                 "--delly",  str(FIX / "delly.vcf"),
                 "--manta",  str(FIX / "manta.vcf"),
                 "--outdir", str(outdir),
                 "--min-informative-pair", "5"],
                capture_output=True, text=True, check=True,
            )
            for f in ("pairwise_relationship_classification.tsv",
                      "pairwise_coefficients.tsv",
                      "candidate_PO_pairs.tsv",
                      "inheritance_segments.tsv",
                      "merged_del_catalogue.json"):
                self.assertTrue((outdir / f).exists(), f"missing {f}")
            # The pairwise table must contain edge_class column.
            header = (outdir / "pairwise_relationship_classification.tsv"
                      ).read_text().splitlines()[0]
            self.assertIn("edge_class", header)
            self.assertIn("theta", header)
            # The candidate PO file must contain at least the two true
            # parent-child pairs (P_F→C and P_M→C). Exclusion may leave
            # direction ambiguous since both directions pass.
            po_lines = (outdir / "candidate_PO_pairs.tsv"
                         ).read_text().splitlines()
            self.assertGreaterEqual(len(po_lines), 1 + 2)  # header + 2 rows
