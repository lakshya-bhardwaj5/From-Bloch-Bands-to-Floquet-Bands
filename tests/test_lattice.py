"""Independent physical benchmarks for the Chapter 2 calculations."""

from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.special import mathieu_a, mathieu_b

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lattice import band_metrics, dispersion, eigensystem, fourier_hopping


class StaticLatticeTests(unittest.TestCase):
    def test_folded_free_particle(self):
        ell = np.arange(-12, 13)
        for k in np.linspace(-1, 1, 17):
            expected = np.sort((k - 2 * ell) ** 2)[:3]
            np.testing.assert_allclose(eigensystem(float(k), 0), expected, atol=2e-12)

    def test_mathieu_edges(self):
        for s in [0, .2, 2, 4, 8, 12]:
            q = s / 4
            center = s / 2 + np.array([mathieu_a(0, q), mathieu_b(2, q), mathieu_a(2, q)])
            edge = s / 2 + np.array([mathieu_b(1, q), mathieu_a(1, q), mathieu_b(3, q)])
            np.testing.assert_allclose(eigensystem(0., s), center, atol=2e-11, rtol=0)
            np.testing.assert_allclose(eigensystem(1., s), edge, atol=2e-11, rtol=0)

    def test_cutoff_convergence_and_inversion(self):
        for s in [0, .25, 2, 4, 8, 12]:
            for k in np.linspace(-1, 1, 21):
                ref = eigensystem(float(k), s, cutoff=16)
                np.testing.assert_allclose(eigensystem(float(k), s), ref, atol=2e-10, rtol=0)
                np.testing.assert_allclose(eigensystem(float(-k), s), ref, atol=2e-10, rtol=0)

    def test_weak_lattice_gap(self):
        s = 1e-3
        self.assertAlmostEqual(band_metrics(s)["edge_gap_ER"] / (s / 2), 1., places=7)

    def test_report_table_2_1(self):
        expected = [[2, .878234, 1.470654, .592420, .996112, .142755],
                    [4, 1.544861, 1.889751, .344890, 1.969357, .085490],
                    [8, 2.486043, 2.609323, .123280, 3.769876, .030799],
                    [12, 3.165608, 3.214620, .049012, 5.304419, .012252]]
        for s, *row in expected:
            met = band_metrics(s)
            values = [met[name] for name in ["E0_center_ER", "E0_edge_ER", "bandwidth_ER", "edge_gap_ER"]]
            values.append(fourier_hopping(s)[1][0])
            np.testing.assert_allclose(values, row, atol=5.1e-7, rtol=0)

    def test_eigenvector_residual(self):
        ell = np.arange(-12, 13)
        for k in [-1., -.37, 0., .44, 1.]:
            energies, vec = eigensystem(k, 4., vectors=True)
            h = np.diag((k - 2 * ell) ** 2 + 2.)
            h += np.diag(np.full(24, -1.), 1) + np.diag(np.full(24, -1.), -1)
            np.testing.assert_allclose(h @ vec, vec * energies, atol=2e-11, rtol=0)
            np.testing.assert_allclose(vec.T @ vec, np.eye(3), atol=2e-12)

    def test_tight_binding_improves(self):
        k = np.linspace(-1, 1, 401)
        errors = []
        for s in [2, 8]:
            mean, hop = fourier_hopping(s)
            energy = dispersion(k, s, bands=1)[:, 0]
            error = np.max(np.abs(energy - mean + 2 * hop[0] * np.cos(np.pi * k))) / np.ptp(energy)
            errors.append(error)
        np.testing.assert_allclose(errors, [.0922, .0107], atol=6e-5, rtol=0)
        self.assertLess(errors[1], errors[0])


if __name__ == "__main__":
    unittest.main()
