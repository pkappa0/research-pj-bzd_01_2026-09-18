import math
import unittest

from src.fetch_gtopdb import normalize_interaction, parse_subunits


class GtoPdbNormalizationTest(unittest.TestCase):
    def test_composition_is_explicit(self):
        parsed = parse_subunits("Human GABA-A receptor alpha1/beta3/gamma2")
        self.assertEqual(parsed["alpha_subunit"], "alpha1")
        self.assertEqual(parsed["beta_subunit"], "beta3")
        self.assertEqual(parsed["gamma_subunit"], "gamma2")
        self.assertEqual(parsed["receptor_composition"], "alpha1/beta3/gamma2")

    def test_missing_fields_are_not_imputed(self):
        row = normalize_interaction(
            {
                "targetId": 123,
                "targetName": "Human GABA-A receptor alpha2/beta3/gamma2",
                "targetType": "LGIC",
                "species": "Homo sapiens",
                "affinityParameter": "Ki",
                "affinityValue": 12.5,
                "affinityUnits": "nM",
                "PMID": "12345678",
            },
            "diazepam",
            "Diazepam",
            "L1",
            "diazepam",
            "2026-01-01T00:00:00+00:00",
            "test",
        )
        self.assertEqual(row["receptor_composition"], "alpha2/beta3/gamma2")
        self.assertEqual(row["affinity_parameter"], "Ki")
        self.assertEqual(row["affinity_value"], 12.5)
        self.assertTrue(math.isnan(row["pKi"]))  # remains missing, not inferred from Ki
        self.assertTrue(row["is_gabaa_related"])


if __name__ == "__main__":
    unittest.main()
