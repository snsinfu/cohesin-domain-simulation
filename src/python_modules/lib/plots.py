import typing

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.inset_locator
import numpy as np
import scipy.ndimage
import scipy.stats

from .kernel_estimates import make_gaussian_kernel, make_epanechnikov_kernel, kernel_regression


def smooth_matrix(matrix: np.ndarray, sigma: float) -> np.ndarray:
    """
    Apply gaussian smoothing on a matrix.
    """
    radius = int(3 * sigma) + 1
    return scipy.ndimage.gaussian_filter(
        matrix, sigma=sigma, radius=radius, mode="reflect",
    )


def kde1d_at(
    xs: np.ndarray,
    **kwargs,
) -> np.ndarray:
    """
    Estimate densities using scipy's gaussian_kde at given points.
    """
    xs = np.asarray(xs)
    valid = np.isfinite(xs)
    data = xs[valid]
    kde = scipy.stats.gaussian_kde(data, **kwargs)
    density = np.full(len(xs), np.nan)
    density[valid] = kde(data)
    return density


def kde2d_at(
    xs: np.ndarray,
    ys: np.ndarray,
    **kwargs,
) -> np.ndarray:
    """
    Estimate densities using scipy's gaussian_kde at given 2D points.
    """
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    valid = np.isfinite(xs) & np.isfinite(ys)
    data = np.array([xs[valid], ys[valid]])
    kde = scipy.stats.gaussian_kde(data, **kwargs)
    density = np.full(len(xs), np.nan)
    density[valid] = kde(data)
    return density


def kde2d(
    x: np.ndarray,
    y: np.ndarray,
    x_at: np.ndarray | None = None,
    y_at: np.ndarray | None = None,
    bandwidth: np.ndarray = np.asarray(1),
    kernel = "gaussian",
) -> np.ndarray:
    """
    Estimate densities using specified kernel at given 2D points.
    """
    match kernel:
        case "gaussian":
            kernel_fn = make_gaussian_kernel(2)
        case "epanechnikov":
            kernel_fn = make_epanechnikov_kernel(2)
        case str(name):
            raise Exception(f"unrecognized kernel {name}")
        case fn:
            kernel_fn = fn

    if x_at is None:
        x_at = x
    if y_at is None:
        y_at = y

    return kernel_regression(
        np.transpose([x_at, y_at]),
        np.transpose([x, y]),
        1,
        kernel_fn=kernel_fn,
        scale=np.asarray(1) / bandwidth,
        average=False,
    )


def update_plot(obj: plt.Axes | plt.Figure):
    """
    Make auto-computed elements (layout, ticks, etc.) in a figure up to date.
    """
    if isinstance(obj, plt.Axes):
        obj = obj.get_figure()
    obj.canvas.draw()


def get_x_grid_line_at(ax: plt.Axes, x: float) -> plt.Line2D | None:
    """
    Get an X-grid Line2D object near specified coordinate.
    """
    for pos, line in _get_x_grid_lines(ax):
        if np.isclose(pos, x):
            return line
    return None


def _get_x_grid_lines(ax: plt.Axes):
    update_plot(ax)
    for line in ax.get_xgridlines():
        match list(line.get_xdata()):
            case [float(x)]:
                yield x, line


def get_y_grid_line_at(ax: plt.Axes, y: float) -> plt.Line2D | None:
    """
    Get a Y-grid Line2D object near specified coordinate.
    """
    for pos, line in _get_y_grid_lines(ax):
        if np.isclose(pos, y):
            return line
    return None


def _get_y_grid_lines(ax: plt.Axes):
    update_plot(ax)
    for line in ax.get_ygridlines():
        match list(line.get_ydata()):
            case [float(y)]:
                yield y, line

def plot_rect(ax: plt.Axes, lx: float, ly: float, tx: float, ty: float, **kwargs) -> mpl.patches.Rectangle:
    """
    Plot a rectangle spanning from lower left (lx, ly) to top right (tx, ty).
    """
    rectangle = mpl.patches.Rectangle((lx, ly), tx - lx, ty - ly, **kwargs)
    ax.add_patch(rectangle)
    return rectangle


def plot_arc(
    ax: plt.Axes,
    cx: float,
    cy: float,
    width: float,
    height: float,
    theta1: float,
    theta2: float,
    angle: float = 0,
    **kwargs,
) ->  mpl.patches.Rectangle:
    """
    Plot an arc of an ellipse centered at (cx, cy) with specified width,
    height, and angle (in degs). The arc starts at angle theta1 and ends at
    angle theta2.
    """
    arc = mpl.patches.Arc(
        (cx, cy),
        width,
        height,
        theta1=theta1,
        theta2=theta2,
        angle=angle,
        **kwargs,
    )
    ax.add_patch(arc)
    return arc


def attach_colorbar(
    ax: plt.Axes,
    sm: mpl.cm.ScalarMappable,
    size: float = 0.03,
    length: float = 0.5,
    margin: float = None,
    orientation: str = "vertical",
    extend: str | None = None,
    outline: bool = True,
):
    """
    Attach a colorbar to a side of an axes.
    """
    if margin is None:
        margin = size

    match orientation:
        case "vertical":
            bbox = (1 + margin, 0, size, length) # right
        case "horizontal":
            bbox = (0.5, -size * 2, 0.5, margin) # bottom
        case _:
            raise ValueError(
                "orientation must be 'vertical' or 'horizontal',"
                + f" not '{orientation}'"
            )

    cbar_ax = mpl_toolkits.axes_grid1.inset_locator.inset_axes(
        ax,
        width="100%", height="100%",
        bbox_to_anchor=bbox, bbox_transform=ax.transAxes,
        borderpad=0,
    )
    fig = ax.get_figure()
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation=orientation, extend=extend)

    if not outline:
        cbar.outline.set_edgecolor("none")

    return cbar


def set_axis(
    ax: plt.Axes,
    spines: str = "none",
    ticks: str | None = None, # default: same as spines
    direction: str | None = None,
):
    """
    Configure axes all at once.
    """
    ALL_SIDES = ["top", "right", "bottom", "left"]

    def map_side_key(key: str) -> list[str]:
        if key == "all":
            return ALL_SIDES
        if key == "none":
            return []
        return [key]

    def make_side_set(sides: str) -> set[str]:
        sides = sides.split()
        sides = set(sum(map(map_side_key, sides), []))
        return sides

    if ticks is None:
        ticks = spines

    tick_options = {}
    if direction is not None:
        tick_options["direction"] = direction

    # Hide everything first.
    for side in ALL_SIDES:
        ax.spines[side].set_visible(False)
        ax.tick_params(which="both", **{side: False, f"label{side}": False})

    # Enable specified sides.
    spines_set = make_side_set(spines)
    for side in spines_set:
        ax.spines[side].set_visible(True)

    ticks_set = make_side_set(ticks)
    for side in ticks_set:
        ax.tick_params(
            which="both",
            **{side: True, f"label{side}": True},
            **tick_options,
        )

    # Move labels when ticks are shown only on one side.
    def specialize(spec, special, default):
        if special in spec and default not in spec:
            return special
        return default

    ax.xaxis.set_label_position(specialize(ticks_set, "top", "bottom"))
    ax.yaxis.set_label_position(specialize(ticks_set, "right", "left"))

    # Ticks should conventionally increase from the labeled side.
    if specialize(ticks_set, "top", "bottom") == "top":
        ax.yaxis.set_inverted(True)
    if specialize(ticks_set, "right", "left") == "right":
        ax.xaxis.set_inverted(True)


def make_scaled_formatter(
    scale: float,
    places: int = 1,
    postfix: str = "",
):
    """
    Create a ticks formatter with constant scaling.
    """
    @mpl.ticker.FuncFormatter
    def formatter(x: float, n: int) -> str:
        return "{0:.{1}f}{2}".format(x * scale, places, postfix)
    return formatter


def register_cmap(cm: mpl.colors.Colormap) -> mpl.colors.Colormap:
    """
    Register colormap to the global namespace.
    """
    if hasattr(plt.cm, "register_cmap"):
        mpl.cm.register_cmap(cmap=cm)
    else:
        mpl.colormaps.register(cm)
    return cm


def lerp_palette(
    nodes: dict[float, object],
    ndiv: int = 20,
) -> np.ndarray:
    """
    Linearly interpolate color nodes to produce a 2D array of RGBA colors.
    """
    points = np.array(sorted(nodes.keys()))
    colors = np.array([mpl.colors.to_rgba(nodes[key]) for key in points])
    x = np.linspace(points[0], points[-1], num=ndiv)
    return np.transpose([np.interp(x, points, colors[:, i]) for i in range(colors.shape[1])])


def faint_color(c: any, level: float = 0.5) -> tuple[float, float, float]:
    """
    Make color more faint.
    """
    r, g, b = plt.matplotlib.colors.to_rgb(c)
    r = (1 - level) * r + level
    g = (1 - level) * g + level
    b = (1 - level) * b + level
    return r, g, b


def reverse_legend(ax: plt.Axes) -> None:
    """
    Reverse the order of legend.
    """
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])


# "fall" colormap with nicer extrema, nice for showing contact heat map.
cm_wrr = register_cmap(
    mpl.colors.LinearSegmentedColormap.from_list(
        "wrr",
        lerp_palette({
            0.00: (1.00, 1.00, 1.00),
            0.01: (1.00, 0.95, 0.90),
            0.50: (1.00, 0.60, 0.00),
            1.00: (0.55, 0.00, 0.00),
        }),
    )
)

# Blue-black-red colormap.
cm_bkr = register_cmap(
    mpl.colors.LinearSegmentedColormap.from_list(
        "bkr", [(0, 0, 1), (0, 0, 0), (1, 0, 0)]
    )
)


# Oklab
class _OkLAB(typing.NamedTuple):
    l: float
    a: float
    b: float


def _oklab(r: float, g: float, b: float) -> _OkLAB:
    # https://bottosson.github.io/posts/oklab/
    l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b

    l_ = l ** (1 / 3)
    m_ = m ** (1 / 3)
    s_ = s ** (1 / 3)

    return _OkLAB(
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_,
    )


def _oklab_distance(x: _OkLAB, y: _OkLAB) -> float:
    lx, ax, bx = x
    ly, ay, by = y
    return ((lx - ly)**2 + (ax - ay)**2 + (bx - by)**2)**(1 / 2)


def _find_nearest_color(
    colors: list[object],
    reference_color: object,
) -> object:
    reference_lab = _oklab(*mpl.colors.to_rgb(reference_color))
    evaluations = [
        (
            color,
            _oklab_distance(_oklab(*mpl.colors.to_rgb(color)), reference_lab),
        )
        for color in colors
    ]
    match min(evaluations, default=None, key=lambda entry: entry[1]):
        case (color, _):
            return color
        case None:
            return None


def _register_cycle_palette():
    # Register named colors that refers to the default color cycle.
    mpl.colors.get_named_colors_mapping().update(
        {
            f"C:{name}": color
            for name, color in _make_cycle_palette(plt.rcParams["axes.prop_cycle"]).items()
        }
    )


def _make_cycle_palette(
    prop_cycle,
    base_colors: list[str] = list("rgbcmykw"),
) -> dict[str, object]:
    prop_colors = [prop["color"] for prop in prop_cycle]
    return {
        name: _find_nearest_color(prop_colors, name) for name in base_colors
    }


_register_cycle_palette()
