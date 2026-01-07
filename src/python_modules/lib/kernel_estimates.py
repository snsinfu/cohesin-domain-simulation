import numpy as np
import scipy.special


def scalar_gaussian_kernel(x):
    return np.sqrt(2 * np.pi) * np.exp(-0.5 * np.square(x))


def scalar_epanechnikov_kernel(x):
    return 0.75 * np.maximum(0, 1 - np.square(x))


def make_gaussian_kernel(dim: int):
    if dim == 1:
        return scalar_gaussian_kernel
    norm = (2 * np.pi) ** (-dim / 2)
    def fn(x):
        return norm * np.exp(-0.5 * np.square(x).sum(-1))
    return fn


def make_epanechnikov_kernel(dim: int):
    if dim == 1:
        return scalar_epanechnikov_kernel
    norm = scipy.special.gamma(dim / 2 + 2) * np.pi**(dim / 2)
    def fn(x):
        return norm * np.maximum(0, 1 - np.square(x).sum(-1))
    return fn


def kernel_regression_at(
    point: np.ndarray,
    sample_points: np.ndarray,
    sample_values: np.ndarray,
    *,
    kernel_fn,
    scale: float | np.ndarray,
    average: bool,
) -> float:
    """
    Kernel regression at given single point.
    """
    weights = kernel_fn((sample_points - point) * scale)
    w_sum = weights.sum(0)
    v_sum = (sample_values * weights).sum(0)
    return v_sum / w_sum if average else v_sum


def kernel_regression(
    points: np.ndarray,
    sample_points: np.ndarray,
    sample_values: np.ndarray,
    *,
    kernel_fn,
    scale: float | np.ndarray = 1,
    average: bool = True,
) -> np.ndarray:
    """
    Kernel regression at given points.
    """
    n = points.shape[0]
    values = np.empty(n)
    for i in range(n):
        values[i] = kernel_regression_at(
            points[i],
            sample_points,
            sample_values,
            kernel_fn=kernel_fn,
            scale=scale,
            average=average,
        )
    return values


def mesh_kernel_regression(
    mesh_spans: tuple,
    sample_points: np.ndarray,
    sample_values: np.ndarray,
    *,
    kernel_fn,
    kernel_coverage: float = 0.05,
    n_div: int = 21,
    average: bool = True,
):
    """
    Kernel regression on a meshgrid.
    """
    scales = np.array([1 / (kernel_coverage * (span[1] - span[0])) for span in mesh_spans])

    coord_points = [np.linspace(*span, num=n_div) for span in mesh_spans]
    grid_coords = np.meshgrid(*coord_points, indexing="xy")

    points = np.reshape(grid_coords, (len(grid_coords), -1)).T
    values = kernel_regression(
        points,
        sample_points,
        sample_values,
        kernel_fn=kernel_fn,
        scale=scales,
        average=average,
    )

    grid_values = values.reshape((n_div, n_div))

    return grid_coords, grid_values


def kernel_flux_estimate(
    points: np.ndarray,
    sequence: np.ndarray,
    *,
    kernel_fn,
    scale: float | np.ndarray = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Estimate flux and density from a sequence of points.

    Args:
        points (N,d): Array of d-dimensional points to estimate flux at.
        sequence (T,d): Sequence of sample points.
        kernel_fn: d-dimensional kernel function to use.
        scale: Scaling factor for normalizing the scale of distance along each
            axis to ~1 (note that the kernel expects dimension-less normalized
            distances). This parameter controls the bandwidth of the kernel.

    Returns:
        A pair (flux, density) of arrays. flux (N,d) gives the estimated flux
        vectors at the given points and density (N) gives the accompanying
        probability density.
    """
    prev_deltas = sequence[None, :-1] - points[:, None]
    next_deltas = sequence[None, 1:] - points[:, None]
    weights = kernel_fn(prev_deltas * scale) * np.prod(scale)
    flux = (weights[:, :, None] * next_deltas).mean(1)
    density = weights.mean(1)
    return flux, density


def mesh_kernel_flux_estimate(
    mesh_spans: tuple,
    sequence: np.ndarray,
    *,
    kernel_fn,
    kernel_coverage: float = 0.05,
    n_div: int = 21,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate flux and density on a meshgrid.
    """
    scale = np.array([
        1 / (kernel_coverage * (span[1] - span[0])) for span in mesh_spans
    ])
    coord_points = [np.linspace(*span, num=n_div) for span in mesh_spans]
    grid_coords = np.meshgrid(*coord_points, indexing="xy")

    points = np.reshape(grid_coords, (len(grid_coords), -1)).T
    flux, density = kernel_flux_estimate(
        points, sequence, kernel_fn=kernel_fn, scale=scale
    )
    grid_flux = flux.reshape((n_div, n_div, len(mesh_spans)))
    grid_density = density.reshape((n_div, n_div))

    return grid_coords, grid_flux, grid_density
