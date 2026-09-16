import numpy as np
import pytest

from inverse_problems.numerics import solve_tridiagonal
from inverse_problems.observations import trapezoidal_integral


def test_solve_tridiagonal_accepts_historical_full_diagonals():
    lower = np.array([0.0, -1.0, -1.0])
    diagonal = np.array([2.0, 2.0, 2.0])
    upper = np.array([-1.0, -1.0, 0.0])
    rhs = np.array([1.0, 0.0, 1.0])

    result = solve_tridiagonal(lower, diagonal, upper, rhs)

    np.testing.assert_allclose(result, [1.0, 1.0, 1.0])


def test_solve_tridiagonal_accepts_compact_diagonals():
    result = solve_tridiagonal([-1.0, -1.0], [2.0, 2.0, 2.0], [-1.0, -1.0], [1.0, 0.0, 1.0])

    np.testing.assert_allclose(result, [1.0, 1.0, 1.0])


def test_solve_tridiagonal_rejects_singular_system():
    with pytest.raises(np.linalg.LinAlgError):
        solve_tridiagonal([], [0.0], [], [1.0])


def test_trapezoidal_integral_matches_linear_function():
    assert trapezoidal_integral([0.0, 1.0, 2.0], 0.5) == pytest.approx(1.0)
