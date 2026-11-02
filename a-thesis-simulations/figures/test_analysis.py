"""Focused checks for diagnostic arithmetic and misleading time interpolation.

Run with: python -m unittest discover -s figures -p 'test_*.py'
"""
import tempfile
import unittest
from pathlib import Path

from publication_style import history_values
from publication_scans import arithmetic, physics_signature, plot_history, read_parameters
import diagnostics as cv
import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg


class DiagnosticTests(unittest.TestCase):
    def test_modal_product_matches_direct_quadrature(self):
        grid = np.array([-2., -.4, .2, 2.])
        a = np.array([[2., -.3], [1.1, .4], [3.2, .7]])
        b = np.array([[1., .1], [1.7, -.2], [2.1, .8]])
        x, w = np.polynomial.legendre.leggauss(3)
        expected = np.sum(np.diff(grid)[:, None] / 2 * w * cv.evaluate(a, x) * cv.evaluate(b, x))
        self.assertAlmostEqual(cv.product_integral(a, b, grid), expected)

    def test_time_mean_is_weighted_by_irregular_spacing(self):
        series = (np.array([0., 1., 10.]), np.array([0., 2., 2.]))
        mean, low, high = cv.window_stats(series, 0., 10.)
        self.assertAlmostEqual(mean, 1.9)
        self.assertEqual((low, high), (0., 2.))

    def test_refuse_extrapolation(self):
        with self.assertRaises(ValueError):
            cv.sample((np.array([0., 1.]), np.ones(2)), [-.1, 1.])

    def test_wall_sum_retains_both_native_grids(self):
        data = {'lo': (np.array([0., .2, 1.]), np.array([0., .4, 2.])),
                'hi': (np.array([0., .7, 1.]), np.array([0., 2.1, 3.]))}
        t, y = history_values(data, [('lo', 1), ('hi', 1)], 0, 1)
        np.testing.assert_array_equal(t, [0., .2, .7, 1.])
        np.testing.assert_allclose(y, 5 * t)

    def test_coarse_moments_do_not_invent_a_15us_transient(self):
        run = dict(uid='baseline', color='blue', linestyle='-',
                   data={'Ei': (np.array([0., 1e-3]), np.array([1., 2.]))})
        fig, ax = plt.subplots()
        try:
            times, values = plot_history(ax, run, [('Ei', 1)], .001 - 15e-6, .001, 1, 1e6)
            np.testing.assert_array_equal(times, [.001])
            np.testing.assert_array_equal(values, [2.])
            self.assertEqual(ax.lines[0].get_linestyle(), 'None')
        finally:
            plt.close(fig)

    def test_sparse_wall_sums_only_plot_shared_samples(self):
        run = dict(uid='baseline', color='blue', linestyle='-', data={
            'lo': (np.array([0., .2, 1.]), np.array([1., 2., 3.])),
            'hi': (np.array([0., .4, 1.]), np.array([4., 5., 6.]))})
        fig, ax = plt.subplots()
        try:
            times, values = plot_history(ax, run, [('lo', 1), ('hi', 1)], .1, .5, 1, 1)
            self.assertEqual(len(times), 0)
            self.assertEqual(len(values), 0)
        finally:
            plt.close(fig)

    def test_phase_parser_keeps_extra_and_continuous_fdp_distinct(self):
        params = dict(alpha_oap='1e-3', alpha_fdp='1', tau_oap='.1', tau_fdp='15e-6',
                      tau_fdp_extra='5*15e-6', num_cycles='5', num_frames_oap='5',
                      num_frames_fdp='5', num_frames_fdp_extra='5*5',
                      cfl_factor_times_omega_max='1/10.0', Z_m='.98')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / 'sim.c').write_text('\n'.join(f'double {k} = {v};' for k,v in params.items()))
            parsed = read_parameters(path)
        self.assertAlmostEqual(parsed['expected_end_s'] - parsed['final_phase_start_s'], 75e-6)
        self.assertAlmostEqual(parsed['expected_end_s'] - parsed['continuous_fdp_start_s'], 90e-6)
        self.assertEqual(parsed['continuous_fdp_start_frame'], 45)
        self.assertEqual(parsed['expected_frame'], 75)

    def test_numeric_input_parser_rejects_code(self):
        self.assertAlmostEqual(arithmetic('1/1000.0'), .001)
        with self.assertRaises(ValueError):
            arithmetic('__import__("os").getcwd()')


class AnnealingTests(unittest.TestCase):
    root = Path(__file__).resolve().parents[1]

    def test_annealing_final_fdp_excludes_the_ramp(self):
        params = read_parameters(self.root / '1x-beams-anneal')
        self.assertEqual(params['expected_frame'], 80)
        self.assertEqual(params['continuous_fdp_start_frame'], 55)
        self.assertAlmostEqual(params['expected_end_s'], .511135)
        self.assertAlmostEqual(params['expected_end_s'] - params['continuous_fdp_start_s'], 75e-6)
        stages = params['annealing_stages']
        self.assertEqual([s['start_frame'] for s in stages], [40, 45, 50, 55])
        np.testing.assert_allclose([s['duration_s'] for s in stages], [.1, .01, .001, 75e-6])

    def test_schedule_review_does_not_hide_physics_changes(self):
        study = 'extra-cycles-fdp-length'
        reference = physics_signature(self.root / '1x-beams', study)
        source = (self.root / '1x-beams-anneal' / 'sim.c').read_text()
        self.assertEqual(physics_signature(self.root / '1x-beams-anneal', study), reference)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for before, after in [('low_pass_damping_rate = 1/5e-6', 'low_pass_damping_rate = 1/50e-6'),
                                  ('{ 1e-3, 1e-2, 1e-1 }', '{ 1e-3, 1e-2, 0.5 }')]:
                self.assertIn(before, source)
                (path / 'sim.c').write_text(source.replace(before, after))
                with self.assertRaisesRegex(ValueError, 'unreviewed annealing phase assembly'):
                    physics_signature(path, study)
            (path / 'sim.c').write_text(source.replace('double Te0 = 940 * eV;', 'double Te0 = 1000 * eV;'))
            self.assertNotEqual(physics_signature(path, study), reference)


class DensityComparisonTests(unittest.TestCase):
    @staticmethod
    def modal(values):
        data = pg.GData(ctx=dict(cells=[2], basis_type='serendipity',
                                poly_order=1, value_form='modal'))
        data.push([np.array([-1., 0., 1.])],
                  pg.gpython.GkylArray.from_numpy(np.asarray(values, dtype=float)))
        return data

    def setUp(self):
        # Uniform computational cells map to physical widths 1.5 and 3.5.
        self.mapping = self.modal([[-2.5 / np.sqrt(2), 1.5 / np.sqrt(6)],
                                   [2.5 / np.sqrt(2), 3.5 / np.sqrt(6)]])
        self.oap = self.modal([[4 * np.sqrt(2), .2], [2 * np.sqrt(2), -.3]])
        self.fdp = self.modal([[np.sqrt(2), -.1], [3 * np.sqrt(2), .4]])

    def test_signed_line_integral_uses_physical_widths(self):
        original = self.oap.get_values().copy()
        result, profile = cv.density_comparison(self.oap, self.fdp, self.mapping)
        # 3*1.5 - 1*3.5 = 1; unweighted computational integration gives 2.
        self.assertAlmostEqual(result['signed_integral_m2'], 1.)
        self.assertAlmostEqual(result['absolute_integral_m2'], 8.)
        self.assertAlmostEqual(result['oap_integral_m2'], 13.)
        self.assertAlmostEqual(result['fdp_integral_m2'], 12.)
        self.assertAlmostEqual(result['relative_to_fdp_percent'], 100 / 12)
        np.testing.assert_allclose(profile[:, 0], [-2., -.5, -.5, 3.])
        np.testing.assert_allclose(profile[:, 3], profile[:, 1] - profile[:, 2])
        self.assertNotEqual(profile[1, 3], profile[2, 3])  # distinct DG traces
        np.testing.assert_array_equal(self.oap.get_values(), original)

    def test_subtraction_sign_and_identical_states(self):
        reverse, _ = cv.density_comparison(self.fdp, self.oap, self.mapping)
        self.assertAlmostEqual(reverse['signed_integral_m2'], -1.)
        self.assertAlmostEqual(reverse['absolute_integral_m2'], 8.)
        same, profile = cv.density_comparison(self.oap, self.oap, self.mapping)
        self.assertEqual(same['signed_integral_m2'], 0.)
        self.assertEqual(same['absolute_integral_m2'], 0.)
        np.testing.assert_array_equal(profile[:, 3], 0.)
        np.testing.assert_array_equal(profile[:, 4], 0.)

    def test_absolute_integral_preserves_zero_mean_cell_variation(self):
        fdp = self.modal([[3 * np.sqrt(2), 0], [3 * np.sqrt(2), 0]])
        oap = self.modal([[3 * np.sqrt(2), np.sqrt(2 / 3)],
                          [3 * np.sqrt(2), 2 * np.sqrt(2 / 3)]])
        result, profile = cv.density_comparison(oap, fdp, self.mapping)
        self.assertAlmostEqual(result['signed_integral_m2'], 0.)
        # int |xi| dxi = 1; physical means are 1/2 and 1.
        self.assertAlmostEqual(result['absolute_integral_m2'], 1.5 / 2 + 3.5)
        np.testing.assert_allclose(profile[:, 4], [-1 / 3, 1 / 3, -2 / 3, 2 / 3])

    def test_absolute_integral_splits_off_center_roots(self):
        fdp = self.modal([[3 * np.sqrt(2), 0], [3 * np.sqrt(2), 0]])
        oap = self.modal([[3.5 * np.sqrt(2), np.sqrt(2 / 3)],
                          [2.5 * np.sqrt(2), np.sqrt(2 / 3)]])
        result, _ = cv.density_comparison(oap, fdp, self.mapping)
        self.assertAlmostEqual(result['absolute_integral_m2'], 5 * .625)

    def test_weak_relative_difference_uses_fdp_denominator(self):
        result, profile = cv.density_comparison(1.2 * self.fdp, self.fdp, self.mapping)
        np.testing.assert_allclose(profile[:, 4], .2, atol=1e-14)
        self.assertAlmostEqual(result['absolute_relative_to_fdp_percent'], 20.)

    def test_nonpositive_fdp_density_is_rejected(self):
        fdp = self.modal([[0, 0], [np.sqrt(2), 0]])
        with self.assertRaisesRegex(ValueError, 'positive FDP density'):
            cv.density_comparison(self.oap, fdp, self.mapping)

    def test_discontinuous_position_mapping_is_rejected(self):
        mapping = self.modal([[-2.5 / np.sqrt(2), 1.5 / np.sqrt(6)],
                              [3.5 / np.sqrt(2), 3.5 / np.sqrt(6)]])
        with self.assertRaisesRegex(ValueError, 'increasing and continuous'):
            cv.density_comparison(self.oap, self.fdp, mapping)


if __name__ == '__main__':
    unittest.main()
