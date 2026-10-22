"""
Nature Journal Figure Style Configuration

Publication-quality figures following Nature journal formatting guidelines
and best practices from Tufte, Rougier et al., and modern data visualization.

Usage:
    from nature_style import set_nature_style, create_figure, save_nature_figure

    set_nature_style()
    fig, ax = create_figure(width='single')
    # ... create plot ...
    save_nature_figure(fig, 'figures/my_figure')

Design Philosophy:
    - Maximize data-ink ratio (Tufte)
    - Avoid chartjunk
    - Design in grayscale first, add color strategically
    - Use direct labeling over legends when possible
    - Maintain typographic hierarchy

References:
    - https://www.nature.com/nature/for-authors/formatting-guide
    - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833
    - https://www.edwardtufte.com/
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# =============================================================================
# Figure Dimensions (Nature specifications)
# =============================================================================

NATURE_WIDTHS = {
    'single': 3.50,      # 89 mm - single column
    '1.5col': 5.04,      # 128 mm - 1.5 columns
    'double': 7.09,      # 180 mm - double column (full width)
}

NATURE_MAX_HEIGHT = 6.69  # 170 mm maximum height


# =============================================================================
# Colorblind-Safe Palettes
# =============================================================================

# Paul Tol's qualitative palette - optimized for colorblind accessibility
# Use for categorical data, limit to 8 colors max for distinguishability
PALETTE_QUALITATIVE = [
    '#332288',  # indigo
    '#88CCEE',  # cyan
    '#44AA99',  # teal
    '#117733',  # green
    '#999933',  # olive
    '#DDCC77',  # sand
    '#CC6677',  # rose
    '#882255',  # wine
    '#AA4499',  # purple
]

# High-contrast two-color schemes (colorblind-safe)
# Choose based on the semantic meaning of your comparison
PALETTE_TWO = {
    'blue_orange': ['#0077BB', '#EE7733'],   # Neutral comparison
    'blue_red': ['#4477AA', '#CC3311'],      # Method vs reference
    'teal_wine': ['#44AA99', '#882255'],     # High contrast
    'blue_gold': ['#2E86AB', '#D4A03C'],     # Elegant contrast
}

# Grayscale for when color isn't needed
PALETTE_GRAY = ['#333333', '#666666', '#999999', '#CCCCCC']

# Recommended colormaps by data type
CMAP_SEQUENTIAL = 'viridis'      # Perceptually uniform, colorblind-safe
CMAP_DIVERGING = 'RdBu_r'        # Symmetric, colorblind-safe
CMAP_SEQUENTIAL_ALT = 'magma'    # Alternative sequential


# =============================================================================
# Style Configuration
# =============================================================================

def set_nature_style(context='paper'):
    """Configure matplotlib for Nature journal style with Tufte-inspired aesthetics.

    Parameters
    ----------
    context : str
        'paper' (default) for publication, 'notebook' for larger display

    Design Principles Applied:
        - Thin lines to maximize data-ink ratio
        - No top/right spines (reduce non-data ink)
        - No gridlines by default
        - Frameless legends
        - Tight layout to minimize whitespace waste
    """
    scale = 1.0 if context == 'paper' else 1.4

    style_dict = {
        # Typography - bumped above Nature's 5-7pt minimum for on-screen/PDF legibility
        # (reviewer feedback: axis labels and tick numbers were too small).
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10 * scale,
        'axes.titlesize': 11 * scale,
        'axes.labelsize': 10 * scale,
        'xtick.labelsize': 9 * scale,
        'ytick.labelsize': 9 * scale,
        'legend.fontsize': 8 * scale,
        'legend.title_fontsize': 9 * scale,
        'figure.titlesize': 11 * scale,

        # Line weights - thin for elegance and data-ink ratio
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'lines.linewidth': 0.75,
        'patch.linewidth': 0.5,
        'hatch.linewidth': 0.5,
        'grid.linewidth': 0.3,

        # Tick settings
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.minor.size': 1.5,
        'ytick.minor.size': 1.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.pad': 2,
        'ytick.major.pad': 2,

        # Tufte-inspired: remove non-data ink
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'axes.axisbelow': True,

        # Legend - minimal, no frame
        'legend.frameon': False,
        'legend.borderpad': 0.3,
        'legend.handlelength': 1.2,
        'legend.handletextpad': 0.4,
        'legend.labelspacing': 0.3,
        'legend.columnspacing': 1.0,

        # Figure settings
        'figure.dpi': 150,
        'savefig.dpi': 450,
        'savefig.format': 'pdf',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.transparent': False,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',

        # Layout - use constrained_layout for consistent spacing
        'figure.constrained_layout.use': True,
        'figure.constrained_layout.h_pad': 0.02,
        'figure.constrained_layout.w_pad': 0.02,
        'figure.autolayout': False,

        # Scatter defaults
        'scatter.edgecolors': 'white',

        # Error bars
        'errorbar.capsize': 2,

        # Math text
        'mathtext.fontset': 'dejavusans',

        # Image interpolation
        'image.interpolation': 'nearest',
        'image.cmap': CMAP_SEQUENTIAL,
    }

    plt.rcParams.update(style_dict)

    # Set default color cycle to colorblind-safe palette
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=PALETTE_QUALITATIVE)

    # Configure seaborn if available
    try:
        import seaborn as sns
        sns.set_palette(PALETTE_QUALITATIVE)
    except ImportError:
        pass


def reset_style():
    """Reset matplotlib to default style."""
    mpl.rcdefaults()


# =============================================================================
# Figure Creation
# =============================================================================

def create_figure(width='single', height=None, nrows=1, ncols=1,
                  aspect_ratio=None, **kwargs):
    """Create a figure with Nature-appropriate dimensions.

    Parameters
    ----------
    width : str or float
        'single' (89mm), '1.5col' (128mm), 'double' (180mm), or inches
    height : float, optional
        Height in inches. If None, calculated from aspect_ratio or golden ratio
    nrows, ncols : int
        Number of subplot rows and columns
    aspect_ratio : float, optional
        Height/width ratio. Defaults to golden ratio (0.618)
    **kwargs
        Additional arguments passed to plt.subplots()

    Returns
    -------
    fig, ax : tuple
        Figure and axes objects
    """
    # Get width in inches
    if isinstance(width, str):
        w = NATURE_WIDTHS.get(width, NATURE_WIDTHS['single'])
    else:
        w = width

    # Calculate height
    if height is None:
        if aspect_ratio is None:
            aspect_ratio = 0.618  # Golden ratio
        # Adjust for multiple rows
        height = w * aspect_ratio * (nrows / ncols) ** 0.5
        height = min(height, NATURE_MAX_HEIGHT)

    return plt.subplots(nrows, ncols, figsize=(w, height), **kwargs)


def create_figure_mm(width_mm, height_mm, nrows=1, ncols=1, **kwargs):
    """Create figure with dimensions specified in millimeters.

    Parameters
    ----------
    width_mm : float
        Width in millimeters
    height_mm : float
        Height in millimeters
    """
    mm_to_inch = 1 / 25.4
    return plt.subplots(nrows, ncols,
                        figsize=(width_mm * mm_to_inch, height_mm * mm_to_inch),
                        **kwargs)


# =============================================================================
# Panel Labels
# =============================================================================

def add_panel_label(ax, label, x=-0.12, y=1.08, **kwargs):
    """Add Nature-style panel label.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    label : str
        The label text (e.g., 'a', 'b', 'c') - use lowercase
    x, y : float
        Position in axes coordinates
    **kwargs
        Additional arguments passed to ax.text()

    Notes
    -----
    Nature style: 8pt, bold, lowercase, upright
    """
    defaults = {
        'fontsize': 9,
        'fontweight': 'bold',
        'va': 'top',
        'ha': 'right',
        'transform': ax.transAxes,
    }
    defaults.update(kwargs)
    ax.text(x, y, label, **defaults)


def add_panel_labels(axes, labels=None, **kwargs):
    """Add panel labels to multiple axes.

    Parameters
    ----------
    axes : array-like of Axes
    labels : list of str, optional
        Custom labels. Defaults to 'a', 'b', 'c', ...
    """
    import string

    if hasattr(axes, 'flatten'):
        axes = axes.flatten()

    if labels is None:
        labels = list(string.ascii_lowercase[:len(axes)])

    for ax, label in zip(axes, labels):
        add_panel_label(ax, label, **kwargs)


# =============================================================================
# Spine and Axis Styling (Tufte-inspired)
# =============================================================================

def despine(ax, top=True, right=True, left=False, bottom=False, trim=False):
    """Remove spines from axes (Tufte-inspired).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    top, right, left, bottom : bool
        Which spines to remove
    trim : bool
        If True, trim remaining spines to data range
    """
    for spine, remove in [('top', top), ('right', right),
                          ('left', left), ('bottom', bottom)]:
        ax.spines[spine].set_visible(not remove)

    if trim:
        # Trim spines to data range
        if not left and ax.get_yticks().size > 0:
            ymin, ymax = ax.get_ylim()
            ax.spines['left'].set_bounds(ymin, ymax)
        if not bottom and ax.get_xticks().size > 0:
            xmin, xmax = ax.get_xlim()
            ax.spines['bottom'].set_bounds(xmin, xmax)


def minimal_spines(ax):
    """Apply minimal spine style - only left and bottom, trimmed."""
    despine(ax, trim=True)


# =============================================================================
# Direct Labeling (reduces legend dependency)
# =============================================================================

def direct_label(ax, x, y, text, offset=(0.02, 0), **kwargs):
    """Add direct label to data point.

    Direct labeling reduces cognitive load compared to legends.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : float
        Data coordinates of the point to label
    text : str
        Label text
    offset : tuple
        (x_offset, y_offset) as fraction of axis range
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_off = offset[0] * (xlim[1] - xlim[0])
    y_off = offset[1] * (ylim[1] - ylim[0])

    defaults = {'fontsize': 6, 'va': 'center', 'ha': 'left'}
    defaults.update(kwargs)

    ax.text(x + x_off, y + y_off, text, **defaults)


def label_line_end(ax, line, text, offset=0.02, **kwargs):
    """Add label at the end of a line (alternative to legend).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    line : matplotlib.lines.Line2D
        The line object to label
    text : str
        Label text
    """
    xdata, ydata = line.get_xdata(), line.get_ydata()
    direct_label(ax, xdata[-1], ydata[-1], f' {text}', offset=(offset, 0), **kwargs)


# =============================================================================
# Statistical Annotations
# =============================================================================

def annotate_correlation(ax, x, y, method='spearman', loc='upper left',
                         show_p=False, **kwargs):
    """Add correlation annotation to scatter plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : array-like
        Data arrays
    method : str
        'spearman' or 'pearson'
    loc : str
        Location: 'upper left', 'upper right', 'lower left', 'lower right'
    show_p : bool
        Whether to show p-value
    """
    from scipy import stats

    if method == 'spearman':
        r, p = stats.spearmanr(x, y)
        symbol = '\u03c1'  # rho
    else:
        r, p = stats.pearsonr(x, y)
        symbol = 'r'

    if show_p:
        if p < 0.001:
            label = f'{method.capitalize()} {symbol} = {r:.2f}\np < 0.001'
        else:
            label = f'{method.capitalize()} {symbol} = {r:.2f}\np = {p:.3f}'
    else:
        label = f'{method.capitalize()} {symbol} = {r:.2f}'

    # Position mapping
    positions = {
        'upper left': (0.05, 0.95, 'top', 'left'),
        'upper right': (0.95, 0.95, 'top', 'right'),
        'lower left': (0.05, 0.05, 'bottom', 'left'),
        'lower right': (0.95, 0.05, 'bottom', 'right'),
    }
    x_pos, y_pos, va, ha = positions.get(loc, positions['upper left'])

    defaults = {'fontsize': 6, 'fontweight': 'bold', 'transform': ax.transAxes,
                'va': va, 'ha': ha}
    defaults.update(kwargs)

    ax.text(x_pos, y_pos, label, **defaults)


def add_significance_bar(ax, x1, x2, y, p_value, height=0.02):
    """Add significance bar between two groups.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x1, x2 : float
        X positions of the two groups
    y : float
        Y position of the bar
    p_value : float
        P-value for annotation
    height : float
        Height of the bar caps as fraction of y-range
    """
    ylim = ax.get_ylim()
    h = height * (ylim[1] - ylim[0])

    # Draw the bar
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y],
            color='black', linewidth=0.5)

    # Add significance text
    if p_value < 0.001:
        text = '***'
    elif p_value < 0.01:
        text = '**'
    elif p_value < 0.05:
        text = '*'
    else:
        text = 'ns'

    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom', fontsize=6)


# =============================================================================
# Reference Lines
# =============================================================================

def add_identity_line(ax, **kwargs):
    """Add identity line (y=x) to axes.

    Useful for method comparison plots.
    """
    defaults = {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5,
                'alpha': 0.7, 'zorder': 0}
    defaults.update(kwargs)

    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
    ]
    ax.plot(lims, lims, **defaults)
    ax.set_xlim(lims)
    ax.set_ylim(lims)


def add_zero_line(ax, orientation='horizontal', **kwargs):
    """Add reference line at zero.

    Parameters
    ----------
    orientation : str
        'horizontal' or 'vertical'
    """
    defaults = {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5,
                'alpha': 0.5, 'zorder': 0}
    defaults.update(kwargs)

    if orientation == 'horizontal':
        ax.axhline(y=0, **defaults)
    else:
        ax.axvline(x=0, **defaults)


# =============================================================================
# Saving Figures
# =============================================================================

def save_nature_figure(fig, filename, formats=None, dpi=450):
    """Save figure in Nature-compatible formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    filename : str
        Base filename without extension
    formats : list of str, optional
        File formats. Defaults to ['pdf', 'png']
    dpi : int
        Resolution for raster formats (Nature minimum: 300, recommended: 450)
    """
    if formats is None:
        formats = ['pdf', 'png']

    for fmt in formats:
        # dpi is passed to all formats: for raster it sets pixel resolution; for
        # PDF/SVG it sets the resolution at which `rasterized=True` artists are
        # rasterized into the vector output. Passing None for PDF causes those
        # artists to be embedded at ~30 dpi as 1-bit indexed PNGs, which loses
        # alpha-blended density (e.g. the dense centre of F5 scatter panels).
        filepath = f'{filename}.{fmt}'
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight',
                    pad_inches=0.02, format=fmt, facecolor='white')
        print(f"Saved: {filepath}")


# =============================================================================
# Color Utilities
# =============================================================================

def get_palette(n_colors=None, palette_type='qualitative'):
    """Get a colorblind-friendly palette.

    Parameters
    ----------
    n_colors : int, optional
        Number of colors (max 9 for qualitative)
    palette_type : str
        'qualitative', 'blue_orange', 'blue_red', 'teal_wine', 'blue_gold', 'gray'

    Returns
    -------
    list : Color hex codes
    """
    palettes = {
        'qualitative': PALETTE_QUALITATIVE,
        'gray': PALETTE_GRAY,
        **PALETTE_TWO
    }

    colors = palettes.get(palette_type, PALETTE_QUALITATIVE)

    if n_colors is not None:
        if n_colors > len(colors):
            raise ValueError(f"Requested {n_colors} colors but palette has {len(colors)}")
        colors = colors[:n_colors]

    return colors


# Backward compatibility alias
get_nature_palette = get_palette


def lighten_color(color, amount=0.3):
    """Lighten a color by mixing with white.

    Parameters
    ----------
    color : str
        Hex color code
    amount : float
        Amount to lighten (0-1)
    """
    import matplotlib.colors as mc
    c = mc.to_rgb(color)
    white = np.array([1, 1, 1])
    return mc.to_hex(np.array(c) * (1 - amount) + white * amount)


def desaturate_color(color, amount=0.3):
    """Desaturate a color (move toward gray).

    Parameters
    ----------
    color : str
        Hex color code
    amount : float
        Amount to desaturate (0-1)
    """
    import matplotlib.colors as mc
    import colorsys

    c = mc.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(*c)
    s = s * (1 - amount)
    return mc.to_hex(colorsys.hls_to_rgb(h, l, s))


# =============================================================================
# Axis Formatting
# =============================================================================

def format_axis_label(text, unit=None):
    """Format axis label following Nature style.

    Nature style: lowercase, first letter capitalized, units in parentheses

    Parameters
    ----------
    text : str
        Label text
    unit : str, optional
        Unit (added in parentheses)

    Returns
    -------
    str : Formatted label

    Examples
    --------
    >>> format_axis_label('concentration', 'mM')
    'Concentration (mM)'
    """
    formatted = text[0].upper() + text[1:].lower() if text else ''
    if unit:
        formatted = f'{formatted} ({unit})'
    return formatted


def set_axis_labels(ax, xlabel=None, ylabel=None, xunit=None, yunit=None):
    """Set axis labels with proper Nature formatting.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    xlabel, ylabel : str
        Label text
    xunit, yunit : str
        Units
    """
    if xlabel:
        ax.set_xlabel(format_axis_label(xlabel, xunit))
    if ylabel:
        ax.set_ylabel(format_axis_label(ylabel, yunit))


# =============================================================================
# Convenience Functions for Common Plot Types
# =============================================================================

def style_scatter(ax, show_identity=False, show_correlation=False, x=None, y=None):
    """Apply common styling to scatter plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    show_identity : bool
        Add identity line (y=x)
    show_correlation : bool
        Show Spearman correlation in corner
    x, y : array-like, optional
        Data for correlation calculation
    """
    if show_identity:
        add_identity_line(ax)

    if show_correlation and x is not None and y is not None:
        annotate_correlation(ax, x, y, method='spearman')


def style_errorbar(ax, orientation='horizontal'):
    """Style for forest plots / error bar plots."""
    add_zero_line(ax, orientation='vertical' if orientation == 'horizontal' else 'horizontal')


# =============================================================================
# Module Demo
# =============================================================================

if __name__ == '__main__':
    # Demo: Create a two-panel figure
    set_nature_style()

    fig, axes = create_figure(width='double', nrows=1, ncols=2)

    # Panel a: scatter plot with correlation
    np.random.seed(42)
    x = np.random.randn(50)
    y = x * 0.8 + np.random.randn(50) * 0.3

    axes[0].scatter(x, y, s=15, alpha=0.7, c=PALETTE_TWO['blue_orange'][0],
                    edgecolors='white', linewidth=0.3)
    add_identity_line(axes[0])
    annotate_correlation(axes[0], x, y, method='spearman')
    set_axis_labels(axes[0], 'method a', 'method b', 'units', 'units')

    # Panel b: bar plot with error bars
    categories = ['A', 'B', 'C', 'D']
    values = [3.2, 4.1, 2.8, 5.0]
    errors = [0.3, 0.4, 0.2, 0.5]
    colors = get_palette(4, 'qualitative')

    axes[1].bar(categories, values, yerr=errors, capsize=2,
                color=colors, edgecolor='white', linewidth=0.5,
                error_kw={'linewidth': 0.5, 'capthick': 0.5})
    set_axis_labels(axes[1], 'category', 'value', yunit='units')
    add_zero_line(axes[1], orientation='horizontal')

    # Add panel labels
    add_panel_labels(axes)

    # Save
    save_nature_figure(fig, '/tmp/claude/nature_demo')
    plt.close()

    print("\nDemo figure saved to /tmp/claude/nature_demo.pdf")
