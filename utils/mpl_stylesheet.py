import matplotlib
from matplotlib import cycler

## Resources from:
## https://matplotlib.org/users/customizing.html
## How to change color and plot styles? 
## https://matplotlib.org/users/dflt_style_changes.html
## matplotlib.rcParams[] = 

def banskt_presentation(black = '#333333', linewidth = 2, ticksize = 8, fontsize = 28, padding = 10, fontfamily = 'latex', colors = 'banskt'):

    if colors == 'banskt':
        mcolors = banskt_colors()
    elif colors == 'kelly':
        mcolors = kelly_colors()

    if fontfamily == 'latex':
        matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[sfdefault,scaled=.85, lining]{FiraSans}', 
                                                      r'\usepackage[cmintegrals]{newtxsf}',
                                                      r'\usepackage{microtype}',
                                                     ]
        matplotlib.rcParams['text.usetex'] = True
    elif fontfamily == 'latex-clearsans':
        matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage[scaled=.86]{ClearSans}',
                                                      r'\usepackage[libertine]{newtxmath}',
                                                      r'\usepackage{microtype}',
                                                     ]
    elif fontfamily == 'system':
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.sans-serif'] = 'DejaVu Sans'
        matplotlib.rcParams['mathtext.fontset'] =  'stixsans'

    # Size
    matplotlib.rcParams['figure.figsize'] = 8, 8
    
    # Fonts
    matplotlib.rcParams['font.size'] = fontsize
    matplotlib.rcParams['text.color'] = black
    matplotlib.rcParams['axes.titlesize'] = fontsize * 1.2
    
    matplotlib.rcParams['axes.labelsize'] = fontsize
    matplotlib.rcParams['axes.labelweight'] = 'normal'
    matplotlib.rcParams['axes.labelcolor'] = black
    
    matplotlib.rcParams['xtick.labelsize'] = fontsize
    matplotlib.rcParams['ytick.labelsize'] = fontsize
    matplotlib.rcParams['legend.fontsize'] = fontsize
    
    # Axes
    matplotlib.rcParams['axes.titlepad'] = 50
    matplotlib.rcParams['axes.edgecolor'] = black
    matplotlib.rcParams['axes.facecolor'] = 'white'
    matplotlib.rcParams['axes.labelpad'] = 20
    matplotlib.rcParams['axes.linewidth'] = linewidth
    
    # Legend
    matplotlib.rcParams['legend.facecolor'] = 'inherit'
    matplotlib.rcParams['legend.edgecolor'] = black
    matplotlib.rcParams['legend.frameon'] = False
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.scatterpoints'] = 1
    matplotlib.rcParams['legend.markerscale'] = 1.0
    # Dimensions as fraction of fontsize
    matplotlib.rcParams['legend.borderpad'] = 0
    matplotlib.rcParams['legend.labelspacing'] = 0.3
    matplotlib.rcParams['legend.handlelength'] = 0.5
    matplotlib.rcParams['legend.handleheight'] = 0.9
    matplotlib.rcParams['legend.handletextpad'] = 0.5
    
    # Ticks
    matplotlib.rcParams['xtick.major.top'] = False
    matplotlib.rcParams['xtick.major.bottom'] = True
    matplotlib.rcParams['xtick.minor.top'] = False
    matplotlib.rcParams['xtick.minor.bottom'] = False
    matplotlib.rcParams['ytick.major.left'] = True
    matplotlib.rcParams['ytick.major.right'] = False
    matplotlib.rcParams['ytick.minor.left'] = False
    matplotlib.rcParams['ytick.minor.right'] = False
    
    matplotlib.rcParams['xtick.major.size'] = ticksize
    matplotlib.rcParams['xtick.minor.size'] = 2 * ticksize / 3.0
    matplotlib.rcParams['ytick.major.size'] = ticksize
    matplotlib.rcParams['ytick.minor.size'] = 2 * ticksize / 3.0
    matplotlib.rcParams['xtick.major.pad'] = padding
    matplotlib.rcParams['xtick.minor.pad'] = padding
    matplotlib.rcParams['ytick.major.pad'] = padding
    matplotlib.rcParams['ytick.minor.pad'] = padding
    matplotlib.rcParams['xtick.major.width'] = linewidth
    matplotlib.rcParams['xtick.minor.width'] = linewidth
    matplotlib.rcParams['ytick.major.width'] = linewidth
    matplotlib.rcParams['ytick.minor.width'] = linewidth
    matplotlib.rcParams['xtick.color'] = black
    matplotlib.rcParams['ytick.color'] = black

    # Color cycle
    matplotlib.rcParams['axes.prop_cycle'] = cycler('color', mcolors)

    # Histogram
    matplotlib.rcParams['hist.bins'] = 20

    # Patches
    # matplotlib.rcParams['patch.facecolor'] = mcolors[0] # doesn't have any effect, comes from prop_cycle
    matplotlib.rcParams['patch.edgecolor'] = black
    matplotlib.rcParams['patch.linewidth'] = linewidth / 2
    matplotlib.rcParams['patch.force_edgecolor'] = True
    
    # For scatter plot, show only left and bottom axes
    matplotlib.rcParams['axes.spines.left'] = True
    matplotlib.rcParams['axes.spines.bottom'] = True
    matplotlib.rcParams['axes.spines.top'] = False
    matplotlib.rcParams['axes.spines.right'] = False

    return

def banskt_colors():
    banskt_colors_hex = [
        '#2D69C4', # blue 
        '#CC2529', # red
        '#93AA00', # Vivid Yellowish Green
        '#535154', # gray
        '#6B4C9A', # purple
        '#FFB300', # Vivid Yellow
        '#922428', # dark brown
        '#948B3D', # olive
        ]
    return banskt_colors_hex

def kelly_colors():
    kelly_colors_hex = [
        '#FFB300', # Vivid Yellow
        '#803E75', # Strong Purple
        '#FF6800', # Vivid Orange
        '#A6BDD7', # Very Light Blue
        '#C10020', # Vivid Red
        '#CEA262', # Grayish Yellow
        '#817066', # Medium Gray

        # The following don't work well for people with defective color vision
        '#007D34', # Vivid Green
        '#F6768E', # Strong Purplish Pink
        '#00538A', # Strong Blue
        '#FF7A5C', # Strong Yellowish Pink
        '#53377A', # Strong Violet
        '#FF8E00', # Vivid Orange Yellow
        '#B32851', # Strong Purplish Red
        '#F4C800', # Vivid Greenish Yellow
        '#7F180D', # Strong Reddish Brown
        '#93AA00', # Vivid Yellowish Green
        '#593315', # Deep Yellowish Brown
        '#F13A13', # Vivid Reddish Orange
        '#232C16', # Dark Olive Green
        ]
    return kelly_colors_hex
