"""Solvers shared by the one-dimensional finite-difference experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _off_diagonal(values: ArrayLike, size: int, name: str, side: str) -> np.ndarray:
    """Normalize an off-diagonal accepted in either common convention."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size not in (size - 1, size):
        raise ValueError(f"{name} must have length n - 1 or n")
    if array.size != size:
        return array
    return array[1:] if side == "lower" else array[:-1]


def solve_tridiagonal(
    lower: ArrayLike,
    diagonal: ArrayLike,
    upper: ArrayLike,
    rhs: ArrayLike,
) -> np.ndarray:
    """Solve a tridiagonal linear system with the Thomas algorithm.

    ``lower`` and ``upper`` may contain either ``n - 1`` active entries or
    ``n`` entries with an unused boundary entry.  The latter convention is
    used by several historical laboratory scripts.
    """
    diag = np.asarray(diagonal, dtype=float)
    vector = np.asarray(rhs, dtype=float)
    if diag.ndim != 1 or vector.ndim != 1 or diag.size != vector.size:
        raise ValueError("diagonal and rhs must be one-dimensional arrays of equal length")
    n = diag.size
    if n == 0:
        raise ValueError("a tridiagonal system cannot be empty")

    lo = _off_diagonal(lower, n, "lower", "lower")
    up = _off_diagonal(upper, n, "upper", "upper")
    modified_diag = diag.copy()
    modified_rhs = vector.copy()

    for index in range(1, n):
        pivot = modified_diag[index - 1]
        if np.isclose(pivot, 0.0):
            raise np.linalg.LinAlgError("zero pivot in tridiagonal system")
        factor = lo[index - 1] / pivot
        modified_diag[index] -= factor * up[index - 1]
        modified_rhs[index] -= factor * modified_rhs[index - 1]

    if np.isclose(modified_diag[-1], 0.0):
        raise np.linalg.LinAlgError("singular tridiagonal system")
    solution = np.empty(n, dtype=float)
    solution[-1] = modified_rhs[-1] / modified_diag[-1]
    for index in range(n - 2, -1, -1):
        if np.isclose(modified_diag[index], 0.0):
            raise np.linalg.LinAlgError("singular tridiagonal system")
        solution[index] = (
            modified_rhs[index] - up[index] * solution[index + 1]
        ) / modified_diag[index]
    return solution
