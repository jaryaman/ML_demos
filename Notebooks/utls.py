import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.ticker
from matplotlib.ticker import FormatStrFormatter

def reset_plots():
    plt.close('all')
    fontsize = 20
    legsize = 15
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('text', usetex=True)
    font = {'size' : fontsize}
    plt.rc('font', **font)
    rc={'axes.labelsize': fontsize,
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'xtick.labelsize':fontsize,
    'ytick.labelsize':fontsize,
    'legend.fontsize': legsize}
    mpl.rcParams.update(**rc)
    mpl.rc('lines', markersize=10)
    plt.rcParams.update({'axes.labelsize': fontsize})
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


def remove_tex_axis(ax, xtick_fmt = '%d', ytick_fmt = '%d'):
	"""
	Makes axes normal font in matplotlib.
	Params:
	xtick_fmt : A string, defining the format of the x-axis
	ytick_fmt : A string, defining the format of the y-axis
	"""
	fmt = matplotlib.ticker.StrMethodFormatter("{x}")
	ax.xaxis.set_major_formatter(fmt)
	ax.yaxis.set_major_formatter(fmt)
	ax.xaxis.set_major_formatter(FormatStrFormatter(xtick_fmt))
	ax.yaxis.set_major_formatter(FormatStrFormatter(ytick_fmt))

def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos.

    pos is an array constructed by packing the meshed arrays of variables
    x_1, x_2, x_3, ..., x_k into its _last_ dimension.

    Source: https://scipython.com/blog/visualizing-the-bivariate-gaussian-distribution/
    """

    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N

def standardize(X):
    """Z-transform an array

    param X: An N x D array where N is the number of examples and D is the number of features

    returns: An N x D array where every column has been rescaled to 0 mean and unit variance
    """
    return (X - X.mean())/X.std(ddof=1)
