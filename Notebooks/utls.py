import matplotlib.pyplot as plt
import matplotlib as mpl

def reset_plots():
	plt.close('all')
	fontsize = 20
	plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rc('text', usetex=True)
	font = {'size' : fontsize}
	plt.rc('font', **font)
	mpl.rc('lines', markersize=10)
	plt.rcParams.update({'axes.labelsize': fontsize})
	mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}', r'\usepackage{amsfonts}']