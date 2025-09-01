import unittest
import numpy as np
from unittest.mock import patch, MagicMock
from safebridge.solvers import NS_Solver


from datetime import timedelta


from src.safebridge.solvers import EW_Solver

class TestNSSolver(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample data that matches the expected structure
        self.sample_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
                'disp': np.array([0.1, 0.15, 0.2, 0.18, 0.12])
            },
            'descending': {
                'ndist': np.array([0.0, 0.3, 0.6, 1.0]),
                'disp': np.array([0.08, 0.14, 0.16, 0.1])
            },
            'deck': {
                'deck_length': np.array([100.0]),
                'span_count': np.array([1])
            }
        }
        
        self.multi_span_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
                'disp': np.array([0.1, 0.12, 0.15, 0.14, 0.11, 0.09])
            },
            'descending': {
                'ndist': np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
                'disp': np.array([0.09, 0.13, 0.16, 0.12, 0.08])
            },
            'deck': {
                'deck_length': np.array([150.0]),
                'span_count': np.array([2])
            }
        }
        
        self.minimal_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.5, 1.0]),
                'disp': np.array([0.1, 0.15, 0.12])
            },
            'descending': {
                'ndist': np.array([0.0, 0.5, 1.0]),
                'disp': np.array([0.08, 0.14, 0.1])
            },
            'deck': {
                'deck_length': np.array([80.0]),
                'span_count': np.array([1])
            }
        }

    def test_init_valid_data(self):
        """Test NS_Solver initialization with valid data."""
        solver = NS_Solver(self.sample_data)
        
        self.assertEqual(solver.data, self.sample_data)
        self.assertIn('ascending', solver.polyfunction)
        self.assertIn('descending', solver.polyfunction)
        self.assertTrue(hasattr(solver.polyfunction['ascending'], '__call__'))
        self.assertTrue(hasattr(solver.polyfunction['descending'], '__call__'))

    def test_init_calls_setup(self):
        """Test that initialization calls setup method."""
        with patch.object(NS_Solver, 'setup') as mock_setup:
            solver = NS_Solver(self.sample_data)
            mock_setup.assert_called_once()

    def test_setup_creates_polyfunction(self):
        """Test that setup method creates polynomial functions correctly."""
        solver = NS_Solver(self.sample_data)
        
        # Check that polyfunction dict is created
        self.assertIsInstance(solver.polyfunction, dict)
        self.assertIn('ascending', solver.polyfunction)
        self.assertIn('descending', solver.polyfunction)
        
        # Check that polynomial functions are numpy poly1d objects
        self.assertIsInstance(solver.polyfunction['ascending'], np.poly1d)
        self.assertIsInstance(solver.polyfunction['descending'], np.poly1d)

    def test_setup_polynomial_degree(self):
        """Test that polynomial functions are degree 2 (quadratic)."""
        solver = NS_Solver(self.sample_data)
        
        # Check polynomial order (degree 2 means order 2)
        self.assertEqual(solver.polyfunction['ascending'].order, 2)
        self.assertEqual(solver.polyfunction['descending'].order, 2)

    def test_quadratic_tilt_ascending(self):
        """Test quadratic_tilt method with ascending data."""
        solver = NS_Solver(self.sample_data)
        tilt = solver.quadratic_tilt('ascending')
        
        self.assertIsInstance(tilt, float)
        self.assertGreaterEqual(tilt, 0)  # Tilt should be non-negative (absolute value)

    def test_quadratic_tilt_descending(self):
        """Test quadratic_tilt method with descending data."""
        solver = NS_Solver(self.sample_data)
        tilt = solver.quadratic_tilt('descending')
        
        self.assertIsInstance(tilt, float)
        self.assertGreaterEqual(tilt, 0)

    def test_quadratic_tilt_calculation(self):
        """Test quadratic_tilt calculation accuracy."""
        solver = NS_Solver(self.sample_data)
        
        # Manual calculation for verification
        xleft = self.sample_data['ascending']['ndist'][0]
        xright = self.sample_data['ascending']['ndist'][-1]
        yleft = solver.polyfunction['ascending'](xleft)
        yright = solver.polyfunction['ascending'](xright)
        expected_tilt = abs(yright - yleft) / self.sample_data['deck']['deck_length'][0]
        
        actual_tilt = solver.quadratic_tilt('ascending')
        self.assertAlmostEqual(actual_tilt, expected_tilt, places=10)

    def test_quadratic_tilt_invalid_keyword(self):
        """Test quadratic_tilt with invalid keyword."""
        solver = NS_Solver(self.sample_data)
        
        with self.assertRaises(KeyError):
            solver.quadratic_tilt('invalid_keyword')

    def test_quadratic_deflection_ascending(self):
        """Test quadratic_deflection method with ascending data."""
        solver = NS_Solver(self.sample_data)
        deflection = solver.quadratic_deflection('ascending')
        
        self.assertIsInstance(deflection, float)
        self.assertGreaterEqual(deflection, 0)

    def test_quadratic_deflection_descending(self):
        """Test quadratic_deflection method with descending data."""
        solver = NS_Solver(self.sample_data)
        deflection = solver.quadratic_deflection('descending')
        
        self.assertIsInstance(deflection, float)
        self.assertGreaterEqual(deflection, 0)

    def test_quadratic_deflection_calculation(self):
        """Test quadratic_deflection calculation accuracy."""
        solver = NS_Solver(self.sample_data)
        
        # Manual calculation for verification
        keyword = 'ascending'
        xleft = self.sample_data[keyword]['ndist'][0]
        xright = self.sample_data[keyword]['ndist'][-1]
        yleft = solver.polyfunction[keyword](xright)
        yright = solver.polyfunction[keyword](xleft)
        
        slope = (yright - yleft) / (xright - xleft)
        intercept = yleft - slope * xleft
        
        xrange = solver._quadratic_x(keyword)
        yinterp = solver.polyfunction[keyword](xrange)
        deflection_vals = yinterp - (slope * xrange + intercept)
        expected_deflection = abs(deflection_vals).max() / self.sample_data['deck']['deck_length'][0]
        
        actual_deflection = solver.quadratic_deflection(keyword)
        self.assertAlmostEqual(actual_deflection, expected_deflection, places=10)

    def test_quadratic_deflection_invalid_keyword(self):
        """Test quadratic_deflection with invalid keyword."""
        solver = NS_Solver(self.sample_data)
        
        with self.assertRaises(KeyError):
            solver.quadratic_deflection('invalid_keyword')

    def test_quadratic_x_ascending(self):
        """Test _quadratic_x method with ascending orbit."""
        solver = NS_Solver(self.sample_data)
        x_vals = solver._quadratic_x('ascending')
        
        self.assertIsInstance(x_vals, np.ndarray)
        self.assertEqual(len(x_vals), 50)  # Default number of points
        self.assertAlmostEqual(x_vals[0], self.sample_data['ascending']['ndist'].min())
        self.assertAlmostEqual(x_vals[-1], self.sample_data['ascending']['ndist'].max())

    def test_quadratic_x_descending(self):
        """Test _quadratic_x method with descending orbit."""
        solver = NS_Solver(self.sample_data)
        x_vals = solver._quadratic_x('descending')
        
        self.assertIsInstance(x_vals, np.ndarray)
        self.assertEqual(len(x_vals), 50)
        self.assertAlmostEqual(x_vals[0], self.sample_data['descending']['ndist'].min())
        self.assertAlmostEqual(x_vals[-1], self.sample_data['descending']['ndist'].max())

    def test_quadratic_x_monotonic(self):
        """Test that _quadratic_x returns monotonically increasing values."""
        solver = NS_Solver(self.sample_data)
        x_vals = solver._quadratic_x('ascending')
        
        # Check that values are monotonically increasing
        self.assertTrue(np.all(np.diff(x_vals) > 0))

    def test_quadratic_y_ascending(self):
        """Test _quadratic_y method with ascending orbit."""
        solver = NS_Solver(self.sample_data)
        y_vals = solver._quadratic_y('ascending')
        
        self.assertIsInstance(y_vals, np.ndarray)
        self.assertEqual(len(y_vals), 50)  # Should match _quadratic_x length

    def test_quadratic_y_descending(self):
        """Test _quadratic_y method with descending orbit."""
        solver = NS_Solver(self.sample_data)
        y_vals = solver._quadratic_y('descending')
        
        self.assertIsInstance(y_vals, np.ndarray)
        self.assertEqual(len(y_vals), 50)

    def test_quadratic_y_consistency(self):
        """Test consistency between _quadratic_x and _quadratic_y."""
        solver = NS_Solver(self.sample_data)
        x_vals = solver._quadratic_x('ascending')
        y_vals = solver._quadratic_y('ascending')
        
        # Manually evaluate polynomial at x_vals
        expected_y = solver.polyfunction['ascending'](x_vals)
        
        np.testing.assert_array_almost_equal(y_vals, expected_y)

    def test_analytical_curve_single_span(self):
        """Test analytical_curve method with single span bridge."""
        solver = NS_Solver(self.sample_data)
        curve = solver.analytical_curve('ascending')
        
        self.assertIsInstance(curve, np.ndarray)
        self.assertGreater(len(curve), 0)

    def test_analytical_curve_multi_span(self):
        """Test analytical_curve method with multi-span bridge."""
        solver = NS_Solver(self.multi_span_data)
        curve = solver.analytical_curve('ascending')
        
        self.assertIsInstance(curve, np.ndarray)
        self.assertGreater(len(curve), 0)

    def test_analytical_curve_insufficient_data(self):
        """Test analytical_curve with insufficient data points."""
        insufficient_data = {
            'ascending': {
                'ndist': np.array([0.0, 1.0]),  # Only 2 points
                'disp': np.array([0.1, 0.12])
            },
            'descending': {
                'ndist': np.array([0.0, 1.0]),
                'disp': np.array([0.08, 0.1])
            },
            'deck': {
                'deck_length': np.array([100.0]),
                'span_count': np.array([1])
            }
        }
        
        solver = NS_Solver(insufficient_data)
        curve = solver.analytical_curve('ascending')
        
        self.assertIsNone(curve)

    def test_analytical_curve_invalid_orbit(self):
        """Test analytical_curve with invalid orbit keyword."""
        solver = NS_Solver(self.sample_data)
        
        with self.assertRaises(KeyError):
            solver.analytical_curve('invalid_orbit')

    @patch('safebridge.solvers.curve_fit')
    def test_analytical_curve_curve_fit_called(self, mock_curve_fit):
        """Test that analytical_curve calls curve_fit for beam displacement."""
        mock_curve_fit.return_value = (np.array([1.0, 2.0, 3.0]), None)
        
        solver = NS_Solver(self.sample_data)
        solver.analytical_curve('ascending')
        
        self.assertTrue(mock_curve_fit.called)

    def test_polynomial_fitting_accuracy(self):
        """Test that polynomial fitting produces reasonable results."""
        # Create data with known quadratic relationship
        x_known = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        y_known = 2 * x_known**2 + 3 * x_known + 1  # Known quadratic: 2x² + 3x + 1
        
        test_data = {
            'ascending': {
                'ndist': x_known,
                'disp': y_known
            },
            'descending': {
                'ndist': x_known,
                'disp': y_known
            },
            'deck': {
                'deck_length': np.array([100.0]),
                'span_count': np.array([1])
            }
        }
        
        solver = NS_Solver(test_data)
        
        # Check that polynomial coefficients are close to expected [2, 3, 1]
        poly_coeffs = solver.polyfunction['ascending'].coefficients
        expected_coeffs = np.array([2.0, 3.0, 1.0])
        
        np.testing.assert_array_almost_equal(poly_coeffs, expected_coeffs, decimal=10)

    def test_data_structure_validation(self):
        """Test behavior with various data structure issues."""
        # Test with missing keys
        incomplete_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.5, 1.0])
                # Missing 'disp' key
            },
            'descending': {
                'ndist': np.array([0.0, 0.5, 1.0]),
                'disp': np.array([0.08, 0.14, 0.1])
            },
            'deck': {
                'deck_length': np.array([100.0]),
                'span_count': np.array([1])
            }
        }
        
        with self.assertRaises(KeyError):
            NS_Solver(incomplete_data)


    def test_minimal_valid_data(self):
        """Test with minimal valid data (3 points for quadratic fit)."""
        solver = NS_Solver(self.minimal_data)
        
        # Should be able to create polynomial functions
        self.assertIsNotNone(solver.polyfunction['ascending'])
        self.assertIsNotNone(solver.polyfunction['descending'])
        
        # Should be able to calculate tilt and deflection
        tilt = solver.quadratic_tilt('ascending')
        deflection = solver.quadratic_deflection('ascending')
        
        self.assertIsInstance(tilt, float)
        self.assertIsInstance(deflection, float)

    def test_large_dataset(self):
        """Test with larger dataset."""
        n_points = 100
        x_large = np.linspace(0, 1, n_points)
        y_large = 0.5 * x_large**2 + 0.1 * x_large + 0.05 + np.random.normal(0, 0.001, n_points)
        
        large_data = {
            'ascending': {
                'ndist': x_large,
                'disp': y_large
            },
            'descending': {
                'ndist': x_large,
                'disp': y_large * 0.9  # Slightly different
            },
            'deck': {
                'deck_length': np.array([200.0]),
                'span_count': np.array([1])
            }
        }
        
        solver = NS_Solver(large_data)
        
        # Should handle large datasets without issues
        tilt = solver.quadratic_tilt('ascending')
        deflection = solver.quadratic_deflection('ascending')
        curve = solver.analytical_curve('ascending')
        
        self.assertIsInstance(tilt, float)
        self.assertIsInstance(deflection, float)
        self.assertIsNotNone(curve)

    def test_numerical_stability(self):
        """Test numerical stability with extreme values."""
        # Test with very small values
        small_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.5, 1.0]),
                'disp': np.array([1e-10, 1.5e-10, 1.2e-10])
            },
            'descending': {
                'ndist': np.array([0.0, 0.5, 1.0]),
                'disp': np.array([0.8e-10, 1.4e-10, 1.0e-10])
            },
            'deck': {
                'deck_length': np.array([100.0]),
                'span_count': np.array([1])
            }
        }
        
        solver = NS_Solver(small_data)
        tilt = solver.quadratic_tilt('ascending')
        deflection = solver.quadratic_deflection('ascending')
        
        self.assertIsInstance(tilt, float)
        self.assertIsInstance(deflection, float)
        self.assertFalse(np.isnan(tilt))
        self.assertFalse(np.isnan(deflection))

    def test_multiple_solver_instances(self):
        """Test creating multiple solver instances."""
        solver1 = NS_Solver(self.sample_data)
        solver2 = NS_Solver(self.multi_span_data)
        
        # Should be independent instances
        self.assertNotEqual(id(solver1.data), id(solver2.data))
        self.assertNotEqual(id(solver1.polyfunction), id(solver2.polyfunction))
        
        # Should produce different results for different data
        tilt1 = solver1.quadratic_tilt('ascending')
        tilt2 = solver2.quadratic_tilt('ascending')
        
        # Results should be different (with high probability)
        self.assertNotAlmostEqual(tilt1, tilt2, places=5)


class TestNSSolverIntegration(unittest.TestCase):
    
    def setUp(self):
        """Set up integration test fixtures."""
        self.integration_data = {
            'ascending': {
                'ndist': np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
                'disp': np.array([0.05, 0.08, 0.12, 0.16, 0.19, 0.20, 0.19, 0.17, 0.14, 0.10, 0.06])
            },
            'descending': {
                'ndist': np.array([0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0]),
                'disp': np.array([0.04, 0.09, 0.14, 0.17, 0.18, 0.16, 0.12, 0.07])
            },
            'deck': {
                'deck_length': np.array([120.0]),
                'span_count': np.array([1])
            }
        }

    def test_complete_workflow(self):
        """Test complete workflow from initialization to results."""
        solver = NS_Solver(self.integration_data)
        
        # Test all major methods work together
        asc_tilt = solver.quadratic_tilt('ascending')
        desc_tilt = solver.quadratic_tilt('descending')
        asc_deflection = solver.quadratic_deflection('ascending')
        desc_deflection = solver.quadratic_deflection('descending')
        asc_curve = solver.analytical_curve('ascending')
        desc_curve = solver.analytical_curve('descending')
        
        # Verify all results are valid
        self.assertIsInstance(asc_tilt, float)
        self.assertIsInstance(desc_tilt, float)
        self.assertIsInstance(asc_deflection, float)
        self.assertIsInstance(desc_deflection, float)
        self.assertIsNotNone(asc_curve)
        self.assertIsNotNone(desc_curve)
        
        # Verify results are reasonable
        self.assertGreater(asc_tilt, 0)
        self.assertGreater(desc_tilt, 0)
        self.assertGreater(asc_deflection, 0)
        self.assertGreater(desc_deflection, 0)

    def test_results_consistency(self):
        """Test that results are consistent across multiple calls."""
        solver = NS_Solver(self.integration_data)
        
        # Call methods multiple times
        tilt1 = solver.quadratic_tilt('ascending')
        tilt2 = solver.quadratic_tilt('ascending')
        deflection1 = solver.quadratic_deflection('ascending')
        deflection2 = solver.quadratic_deflection('ascending')
        
        # Results should be identical
        self.assertEqual(tilt1, tilt2)
        self.assertEqual(deflection1, deflection2)

    def test_cross_orbit_comparison(self):
        """Test comparing results between ascending and descending orbits."""
        solver = NS_Solver(self.integration_data)
        
        asc_tilt = solver.quadratic_tilt('ascending')
        desc_tilt = solver.quadratic_tilt('descending')
        asc_deflection = solver.quadratic_deflection('ascending')
        desc_deflection = solver.quadratic_deflection('descending')
        
        # Results should be in reasonable range relative to each other
        self.assertLess(abs(asc_tilt - desc_tilt) / max(asc_tilt, desc_tilt), 2.0)  # Within 200%
        self.assertLess(abs(asc_deflection - desc_deflection) / max(asc_deflection, desc_deflection), 2.0)

class TestEWSolver(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample time overlap information
        start_date = np.datetime64('2020-01-01')
        end_date = np.datetime64('2020-03-01')
        
        asc_dates = np.array([
            start_date - timedelta(days=10),
            start_date - timedelta(days=5),
            start_date + timedelta(days=0),
            start_date + timedelta(days=15),
            start_date + timedelta(days=30),
            start_date + timedelta(days=45),
            end_date + timedelta(days=10)
        ], dtype='datetime64[D]')
        
        dsc_dates = np.array([
            start_date - timedelta(days=8),
            start_date - timedelta(days=2),
            start_date + timedelta(days=5),
            start_date + timedelta(days=20),
            start_date + timedelta(days=40),
            start_date + timedelta(days=60),
            end_date + timedelta(days=15)
        ], dtype='datetime64[D]')
        
        self.time_overlap_info = {
            'rmin': start_date,
            'rmax': end_date,
            'ascending': {'date': asc_dates},
            'descending': {'date': dsc_dates}
        }
        
        # Test parameters
        self.theta_asc = 35.0  # incidence angle
        self.theta_dsc = 38.0  # incidence angle
        self.alpha_asc = 100.0  # azimuth angle
        self.alpha_dsc = 260.0  # azimuth angle
        
        # Initialize solver
        self.solver = EW_Solver(
            self.time_overlap_info,
            self.theta_asc,
            self.theta_dsc,
            self.alpha_asc,
            self.alpha_dsc
        )
        
        # Sample displacement data
        self.ascending_ts = np.array([0.0, 1.0, 2.0, 4.0, 7.0, 9.0, 12.0])
        self.descending_ts = np.array([0.5, 1.5, 3.0, 5.0, 8.0, 10.0, 13.0])
        
        # Sample bridge parameters
        self.bridge_azimuth = 45.0
        self.deck_length = 100.0
        
    def test_initialization(self):
        """Test if EW_Solver initializes correctly."""
        self.assertIsInstance(self.solver, EW_Solver)
        self.assertEqual(self.solver.timeOverlapInfo, self.time_overlap_info)
        self.assertAlmostEqual(self.solver.theta_asc, np.deg2rad(self.theta_asc))
        self.assertAlmostEqual(self.solver.theta_dsc, np.deg2rad(self.theta_dsc))
        self.assertAlmostEqual(self.solver.alpha_asc, np.deg2rad(self.alpha_asc - 90))
        self.assertAlmostEqual(self.solver.alpha_dsc, np.deg2rad(self.alpha_dsc - 90))
        
    def test_process_time_overlap_info(self):
        """Test if time overlap info is processed correctly."""
        # Check masks
        start_date = self.time_overlap_info['rmin']
        end_date = self.time_overlap_info['rmax']
        asc_time = self.time_overlap_info['ascending']['date']
        dsc_time = self.time_overlap_info['descending']['date']
        
        asc_mask = (asc_time >= start_date) & (asc_time <= end_date)
        dsc_mask = (dsc_time >= start_date) & (dsc_time <= end_date)
        
        # Access private attributes via mangled names
        self.assertTrue(np.array_equal(self.solver._EW_Solver__asc_mask, asc_mask))
        self.assertTrue(np.array_equal(self.solver._EW_Solver__dsc_mask, dsc_mask))
        
    def test_ts_interpolation(self):
        """Test time series interpolation method."""
        # Extract limited data for testing
        asc_dates = self.time_overlap_info['ascending']['date']
        asc_mask = (asc_dates >= self.time_overlap_info['rmin']) & (asc_dates <= self.time_overlap_info['rmax'])
        asc_num = (asc_dates[asc_mask] - self.time_overlap_info['rmin']).astype('timedelta64[D]').astype(int)
        
        # Take a subset of the time series data
        test_ts = np.array([2.0, 4.0, 7.0, 2.0, 4.0, 7.0, 1.0])  # Values at days 0, 15, 30
        
        # Call ts_interpolation
        interpolated = self.solver.ts_interpolation(test_ts, asc_num)
        
        # Verify it's a numpy array with expected size
        self.assertIsInstance(interpolated, np.ndarray)
        
        # Verify some values are interpolated correctly
        # Day 0 should have value 2.0
        self.assertAlmostEqual(interpolated[0], 2.0, places=5)
        
        # Find the index for day 15 in the combined dates
        combined_dates = self.solver._EW_Solver__combined_dates
        day15_idx = np.where(combined_dates == 15)[0][0]
        # Day 15 should have value 4.0
        self.assertAlmostEqual(interpolated[day15_idx], 4.0, places=5)

        
    def test_los_long_vert_displacement(self):
        """Test LOS to longitudinal and vertical displacement calculation."""
        # Create sample interpolated displacements
        interp_asc_disp = np.array([0.0, 1.0, 2.0, 3.0])
        interp_dsc_disp = np.array([0.0, 1.5, 3.0, 4.5])
        
        # Call the method
        dL, dV = self.solver.los_long_vert_displacement(
            interp_asc_disp, 
            interp_dsc_disp, 
            self.bridge_azimuth
        )
        
        # Check output types and shapes
        self.assertIsInstance(dL, np.ndarray)
        self.assertIsInstance(dV, np.ndarray)
        self.assertEqual(dL.shape, interp_asc_disp.shape)
        self.assertEqual(dV.shape, interp_dsc_disp.shape)
        
        # Verify the calculation works - we don't test exact values since they depend on 
        # complex trigonometric calculations, but we can check they're not all zeros or NaNs
        self.assertFalse(np.allclose(dL, 0.0))
        self.assertFalse(np.allclose(dV, 0.0))
        self.assertFalse(np.any(np.isnan(dL)))
        self.assertFalse(np.any(np.isnan(dV)))
        
    def test_get_tilt(self):
        """Test get_tilt calculation."""
        # Test case 1: Valid data with N and S keys
        data_store = {
            'N': {'vert': 5.0},
            'S': {'vert': 2.0}
        }
        tilt = self.solver.get_tilt(data_store, self.deck_length)
        self.assertAlmostEqual(tilt, 3.0 / self.deck_length)
        
        # Test case 2: Missing keys
        incomplete_data = {
            'N': {'vert': 5.0},
            'C': {'vert': 3.0}
        }
        tilt = self.solver.get_tilt(incomplete_data, self.deck_length)
        self.assertIsNone(tilt)
        
        # Test case 3: NaN values
        nan_data = {
            'N': {'vert': np.nan},
            'S': {'vert': 2.0}
        }
        tilt = self.solver.get_tilt(nan_data, self.deck_length)
        self.assertIsNone(tilt)
        
    def test_get_deflection(self):
        """Test get_deflection calculation."""
        # Test case 1: Valid data with N, C, and S keys
        ndist = np.array([0.0, 0.5, 1.0])
        data_store = {
            'N': {'long': 2.0},
            'C': {'long': 5.0},
            'S': {'long': 3.0}
        }
        deflection = self.solver.get_deflection(data_store, ndist, self.deck_length)
        self.assertIsNotNone(deflection)
        self.assertGreater(deflection, 0)
        
        # Test case 2: Missing keys
        incomplete_data = {
            'N': {'long': 2.0},
            'S': {'long': 3.0}
        }
        deflection = self.solver.get_deflection(incomplete_data, ndist, self.deck_length)
        self.assertIsNone(deflection)
        
        # Test case 3: NaN values
        nan_data = {
            'N': {'long': np.nan},
            'C': {'long': 5.0},
            'S': {'long': 3.0}
        }
        deflection = self.solver.get_deflection(nan_data, ndist, self.deck_length)
        self.assertIsNone(deflection)


if __name__ == '__main__':
    unittest.main()