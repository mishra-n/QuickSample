import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import Table
from astropy.io import ascii
import matplotlib.image as mpimg
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord, Distance
from astropy import units as u
from astropy.constants import c, M_sun, G

#import redshifting as redshift
import astropy
from astropy.io import fits
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.nddata import NDData
from astropy.constants import m_e, c, e
from astropy.nddata import StdDevUncertainty
from astropy.convolution import Gaussian1DKernel, Box1DKernel, convolve

import scipy.interpolate as interp
from scipy.integrate import simps
from scipy.ndimage import gaussian_filter1d, uniform_filter1d
import scipy.stats as stats
import matplotlib as mpl
from astropy.table import QTable, Table, vstack
import corner
import mpdaf
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_lines
from specutils.fitting import fit_generic_continuum

from mendeleev import element
import pandas as pd
import time

mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams['figure.facecolor'] = 'white'
plt.rc('text', usetex=False)
plt.rc('font', family='serif',size=20)
plt.rc('axes', linewidth=1.5)
plt.rc('axes', labelsize=25)
plt.rc('xtick', labelsize=15, direction='in')
plt.rc('ytick', labelsize=15, direction='in')
plt.rc('xtick', top=True)
plt.rc('ytick', right=True)
plt.rc('xtick.minor', visible=True)
plt.rc('ytick.minor', visible=True)
plt.rc('xtick.major',size=10, pad=4)
plt.rc('xtick.minor',size=5, pad=4)
plt.rc('ytick.major',size=10)
plt.rc('ytick.minor',size=5)
plt.rc('legend', fontsize=15)
