"""Small, reusable observation operators."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def trapezoidal_integral(values: ArrayLike, step: float) -> float:
    """Return the composite trapezoidal integral on a uniform grid."""
    samples = np.asarray(values, dtype=float)
    if samples.ndim != 1 or samples.size < 2:
        raise ValueError("at least two one-dimensional samples are required")
    if step <= 0:
        raise ValueError("step must be positive")
    return float(step * (samples.sum() - 0.5 * (samples[0] + samples[-1])))
