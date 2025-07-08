import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import rc, rcParams
#from matplotlib.font_manager import fontManager, FontProperties
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
#rc('font',**{'family':'serif','serif':['Computer Modern']})
rcParams.update({'font.size': 14})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = 'on'
rcParams['ytick.right'] = 'on'
rcParams['axes.xmargin'] = 0
rcParams['xtick.major.size'] = 4
rcParams['xtick.major.width'] = 1.2

def readspecfits(fitsfile):
    # Read a spectrum
    hdu0 = fits.open(fitsfile)[0]
    
    freq_spw = (np.arange(hdu0.header['NAXIS3'])+1.0-hdu0.header['CRPIX3' ]) * hdu0.header['CDELT3'] + hdu0.header['CRVAL3']
    freq_spw = freq_spw * 1e-9 # GHz
 
    spec_spw = hdu0.data[:,0,0]
    bmaj, bmin = hdu0.header['BMAJ']*3600, hdu0.header['BMIN']*3600 # arcsec
    freq = hdu0.header['RESTFRQ'] * 1e-9 # GHz
    # Convert intensty to Tb
    spec_spw = spec_spw * 1.224e6 / (freq**2) / (bmaj * bmin)

    return freq_spw, spec_spw

# Read spectra
freq_spw33, spec_spw33 = readspecfits('SgrC_coreMM1_average_uid___A001_X15a0_X172.s38_0.Sgr_A_star_sci.spw33.cube.I.iter1.image.pbcor.fits')
freq_spw35, spec_spw35 = readspecfits('SgrC_coreMM1_average_uid___A001_X15a0_X172.s38_0.Sgr_A_star_sci.spw35.cube.I.iter1.image.pbcor.fits')

freq_spw33_brick, spec_spw33_brick = readspecfits('BrickMaserCore_average_uid___A001_X15a0_X190.s38_0.Sgr_A_star_sci.spw33.cube.I.iter1.image.pbcor.fits')
spec_spw33_brick = np.interp(freq_spw33, freq_spw33_brick, spec_spw33_brick)

freq_spw35_brick, spec_spw35_brick = readspecfits('BrickMaserCore_average_uid___A001_X15a0_X190.s38_0.Sgr_A_star_sci.spw35.cube.I.iter1.reclean.image.pbcor.fits')
spec_spw35_brick = np.interp(freq_spw35, freq_spw35_brick, spec_spw35_brick)

# Plot spectra
fig = plt.figure(221,figsize=(14, 8))

fig.text(0.005,0.5,r'$T_B$ (K)',rotation='vertical')
fig.text(0.45,0.01,r'Sky Frequency (GHz)')

ylim0, ylim1 = -0.9, 16.1

# The first spectrum
spec1 = fig.add_subplot(221,position=[0.04,0.53,0.465,0.46])
spec1.step(freq_spw33,spec_spw33,color='blue',linewidth=0.7)
spec1.set_ylim(ylim0, ylim1)
spec1.spines['right'].set_visible(False)
spec1.yaxis.set_ticks_position('left') 
spec1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
spec1.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
spec1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
spec1.yaxis.set_major_locator(ticker.MultipleLocator(5))
spec1.set_xticklabels(())
spec1.grid(True, alpha=0.2, ls=':', color='k')
spec1.text(97.75,ylim1-2.0,r'Sgr C MM1 hot core ($V_\mathrm{lsr}$=$-$50 km s$^{-1}$)',color='black',size='large',weight='black')

# The second spectrum
spec2 = fig.add_subplot(222,position=[0.53,0.53,0.465,0.46])
spec2.step(freq_spw35,spec_spw35,color='blue',linewidth=0.7)
spec2.set_ylim(ylim0, ylim1)
#spec2.yaxis.set_visible(False)
spec2.spines['left'].set_visible(False)
#spec2.spines['right'].set_visible(False)
spec2.yaxis.set_ticks_position('right')
spec2.yaxis.set_ticklabels([])
spec2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
spec2.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
spec2.yaxis.set_minor_locator(ticker.MultipleLocator(1))
spec2.yaxis.set_major_locator(ticker.MultipleLocator(5))
spec2.grid(True, alpha=0.2, ls=':', color='k')
spec2.set_xticklabels(())

# The third spectrum
spec3 = fig.add_subplot(223,position=[0.04,0.07,0.465,0.46])
spec3.step(freq_spw33,spec_spw33_brick,color='blue',linewidth=0.7)
spec3.set_ylim(ylim0, ylim1)
spec3.spines['right'].set_visible(False)
spec3.yaxis.set_ticks_position('left') 
spec3.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
spec3.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
spec3.yaxis.set_minor_locator(ticker.MultipleLocator(1))
spec3.yaxis.set_major_locator(ticker.MultipleLocator(5))
spec3.grid(True, alpha=0.2, ls=':', color='k')
#spec3.set_xticklabels(())
spec3.text(97.75,ylim1-2.0,r'G0.253+0.016 H$_2$O-maser core ($V_\mathrm{lsr}$=40 km s$^{-1}$)',color='black',size='large',weight='black')

# The fourth spectrum
spec4 = fig.add_subplot(224,position=[0.53,0.07,0.465,0.46])
spec4.step(freq_spw35,spec_spw35_brick,color='blue',linewidth=0.7)
spec4.set_ylim(ylim0, ylim1)
#spec4.yaxis.set_visible(False)
spec4.spines['left'].set_visible(False)
#spec4.spines['right'].set_visible(False)
spec4.yaxis.set_ticks_position('right')
spec4.yaxis.set_ticklabels([])
spec4.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
spec4.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
spec4.yaxis.set_minor_locator(ticker.MultipleLocator(1))
spec4.yaxis.set_major_locator(ticker.MultipleLocator(5))
spec4.grid(True, alpha=0.2, ls=':', color='k')
#spec4.set_xticklabels(())


# Dictionaries for the lines
lines_lsb = {
    r'CS 2-1': 97.98095,
    r'$^{34}$SO 3$_2$-2$_1$': 97.7154,
    r'CH$_3$CHO 5$_{1,4}$-4$_{1,3}$ A': 98.90095,
    r'H40$\alpha$': 99.02295,
    r'SO 3$_2$-2$_1$': 99.29987,
    r'CH$_3$CHO 5$_{1,4}$-4$_{1,3}$ E': 98.863314,
}

lines_usb = {
    r'HC$_3$N 11-10': 100.0763,
}

##### Hot core spectra labels
vlsr = -50 # km/s
for line,freq in lines_lsb.items():
    freq *= (1 - vlsr/3e5)
    spec1.axvline(x=freq,linestyle='dashed',color='r',linewidth=0.6)
    if line == 'CH$_3$CHO 5$_{1,4}$-4$_{1,3}$ A':
        spec1.text(freq+0.015,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)
    else:
        spec1.text(freq+0.005,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)

for line,freq in lines_usb.items():
    freq *= (1 - vlsr/3e5)
    spec2.axvline(x=freq,linestyle='dashed',color='r',linewidth=0.6)
    spec2.text(freq+0.005,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)


##### Brick core spectra labels

vlsr = 40 # km/s
for line,freq in lines_lsb.items():
    freq *= (1 - vlsr/3e5)
    spec3.axvline(x=freq,linestyle='dashed',color='r',linewidth=0.6)
    if line == 'CH$_3$CHO 5$_{1,4}$-4$_{1,3}$ A':
        spec3.text(freq+0.015,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)
    else:
        spec3.text(freq+0.005,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)

for line,freq in lines_usb.items():
    freq *= (1 - vlsr/3e5)
    spec4.axvline(x=freq,linestyle='dashed',color='r',linewidth=0.6)
    spec4.text(freq+0.005,4,line,rotation=90,color='black',size='medium',weight=1000,ha='center',va='bottom',zorder=90)


####
fig.savefig('aces_spw3335_spec_useAdamspectra.pdf',dpi=300)
