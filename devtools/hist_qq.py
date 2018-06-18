import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import chi2
from scipy.stats import uniform
from scipy import special

def normal_quantile_plot(data, mu = 0, sigma = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = norm.ppf(mq, loc = mu, scale = sigma)
    # x = mu + sigma * np.sqrt(2) * special.erfinv(2 * mq - 1) # qunatile function for normal distribution (see wiki)
    return x, y

def normal_qq(data, mu, sigma):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = stats.norm.ppf(mq, loc = mu, scale = sigma)
    return x, y

def uniform_quantile_plot_old(data, a = 0, b = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = a + mq * (b - a) # quantile function for uniform distribution: a + p(b - a)
    return x, y

def uniform_quantile_plot(data, a = 0, b = 1):
    y = data[np.argsort(data)]
    xrand = np.random.uniform(0, 1, size=y.shape[0])
    x = xrand[np.argsort(xrand)]
    return x, y


def chisquare_quantile_plot(data, df = 1):
    y = data[np.argsort(data)]
    mq = (np.arange(y.shape[0] + 2) / (y.shape[0] + 1))[1:-1]
    x = chi2.ppf(mq, df) # quantile function for chisquare distribution from python
    return x, y

def get_pdf(dist, x, df = 1, loc = 0, scale = 1):
    if dist == 'normal':
        y = norm.pdf(x, loc = loc, scale = scale)
    elif dist == 'uniform':
        y = uniform.pdf(x, loc = loc, scale = scale)
    elif dist == 'chi2':
        y = chi2.pdf(x, df)
    else:
        print ("No distribution found with name {:s}".format(dist))
        y = np.zeros_like(x)
    return y

def get_quantile(dist, data, df = 1, loc = 0, scale = 1):
    if dist == 'normal':
        x, y = normal_quantile_plot(data, mu = loc, sigma = scale)
    elif dist == 'uniform':
        x, y = uniform_quantile_plot(data)
        x = - np.log10(x)
        y = - np.log10(y)
    elif dist == 'chi2':
        x, y = chisquare_quantile_plot(data, df = df)
    else:
        print ("No distribution found with name {:s}".format(dist))
        y = np.zeros_like(x)
    return x, y

def plot(ax1, ax2, data, nbins, dist, df = 0, loc = 0, scale = 1, size = 1, label = 'Simulation', xlim = None, trans = None):

    axisfontsize = 20 * size
    labelfontsize = 15 * size
    padwidth = 10 * size
    msize = 10 * size
    wsize = 2 * size
    ticklen = 5 * size
    borderwidth = 2 * size
    bordercolor = 'black'
    nsnps = data.shape[0]

    banskt_colors_hex = [
        '#2D69C4', # blue 
        '#FFB300', # Vivid Yellow
        '#93AA00', # Vivid Yellowish Green
        '#CC2529', # red
        '#535154', # gray
        '#6B4C9A', # purple
        '#922428', # dark brown
        '#948B3D', # olive
        ]
    colors = banskt_colors_hex

    if trans is None:
        trans = np.zeros(nsnps).astype(bool)

    xmin = 0 if dist == 'uniform' else np.min(data)
    xmax = 1 if dist == 'uniform' else np.max(data)
    bins = np.linspace(xmin, xmax, nbins)
    xvals = [(bins[i] + bins[i+1]) / 2 for i in range(nbins - 1)]
    h, _b = np.histogram(data, bins=bins, density=False)
    _f = np.sum(h * np.diff(_b))

    h1, _ = np.histogram(data[~trans], bins=bins, density=False)
    h1 = h1 / _f
    ax1.fill_between(xvals, h1, 0, color=colors[7], alpha = 0.2, label = label)
    if sum(trans) > 0:
        h2, _ = np.histogram(data[trans], bins=bins, density=False)
        h2 = h2 / _f
        ax1.fill_between(xvals, h2, 0, color=colors[3], alpha = 0.7, label = label)
    
    xvals = np.linspace(xmin, xmax, data.shape[0])
    yvals = get_pdf(dist, xvals, df = df, loc = loc, scale = scale)
    ax1.plot(xvals, yvals, lw = wsize, color=colors[7], label = 'Theory')

###     legend_elements = [Line2D([0], [0], color='b', lw=4, label='Line'),
###                    Line2D([0], [0], marker='o', color='w', label='Scatter',
###                           markerfacecolor='g', markersize=15),
###                    Patch(facecolor='orange', edgecolor='r',
###                          label='Color Patch')]
###     ax.legend(handles=legend_elements, loc='center')
##    
###     mh, ml = ax1.get_legend_handles_labels()
##    legend = ax1.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98),
##                        handlelength = 3.0,
##                        handletextpad = 1.0,
##                        markerscale=5,
##                        ncol = 1,
##                        frameon = True, borderpad = 1.5, labelspacing = 1.5
##                        #title = legendtitle
##                       )
###     for l in legend.legendHandles:
###         l.set_alpha(1)
##
##    lframe = legend.get_frame()
##    lframe.set_edgecolor(bordercolor)
##    lframe.set_linewidth(borderwidth)
##    for fonts in ([legend.get_title()] + legend.texts):
##        fonts.set_fontsize(labelfontsize)
##        fonts.set_color(bordercolor)

    x, y = get_quantile(dist, data, df = df, loc = loc, scale = scale)
    xmin = min(np.min(x), np.min(y))
    xmax = max(np.max(x), np.max(y))
    ax2.scatter(x, y, s = msize, color=colors[3], alpha = 0.5)
    ax2.plot([xmin, xmax], [xmin, xmax], lw = wsize / 4, ls = 'dashed', color=colors[4])

    ax1.set_xlabel('x', {'size': axisfontsize}, labelpad = padwidth)
    ax1.set_ylabel('PDF', {'size': axisfontsize}, labelpad = padwidth)
    
    ax2.set_xlabel('Expected', {'size': axisfontsize}, labelpad = padwidth)
    ax2.set_ylabel('Observed', {'size': axisfontsize}, labelpad = padwidth)

    xticks = ax2.get_xticks()
    ax2.set_yticks(xticks)
    ax2.set_xticks(xticks)
    if xlim is not None:
        ax1.set_xlim(xlim)
    for ax in [ax1, ax2]:
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        ax.tick_params(axis='both', which = 'major',
                       length = ticklen, width = borderwidth, pad=padwidth,
                       labelsize = labelfontsize,
                       color = bordercolor,
                       labelcolor = bordercolor,
                       bottom = True, top = False, left = True, right = False,
                      )
        for side, border in ax.spines.items():
            border.set_linewidth(borderwidth)
            border.set_color(bordercolor)
            
    font_properties = {'family':'sans-serif', 'weight': 'bold', 'size': labelfontsize}
    
    return None
