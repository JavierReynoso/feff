#Code to compute the efficiency function for a certain energy injecting dark matter particle model

"""
    
Created on April 23 19:00  2017
@author J. Reynoso-Cordova

eventhought this code is of the author's property you must cite 1506.03811
this is only to make a Python version.

The necesarry files to run this code are in

https://faun.rc.fas.harvard.edu/epsilon/detaileddeposition/annihilation/

you must either have the files and the code on the same Directory or change it in section II.

Several spectrums (arXiv:1705.00777v1)

"""

from pylab import*
from matplotlib import*
from math import*
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from astropy.io import fits

#############################################################
##    Define the electron-positron or photon spectrum      ##
#############################################################

AssumedIonizationHistory ="3keV"        #you can change to "SSCK", "3keV" to use planck constraints

"""
    
    Here you either import the elec_pos/photon dN/dE and build an interpolation function or define the function
    
    
    """

#Electron-Positron Spectrum for Dark matter annihilating into Charged Pions

def Energy_Spectrum_electrons(eng,mass):
        return 0.

masspion=135.e6
def logeng_min(mass):
    return log(masspion**2/(4.*mass))

def logeng_max(mass):
    return mass

def Energy_Spectrum_photons(eng,mass):
    return 4./(mass - masspion**2/(4.*mass))



############################################################
##      Import the f_eff^{e^+} and f_eff^{\gamma}         ##
############################################################


"""
    There is nothing else to change here except the directory in the open file
"""


##     feff for photons   ##


Species_1="phot"
data_feff_fits_phot=fits.open('feff_summary_'+Species_1+'_'+AssumedIonizationHistory+'.fits')
Grid_phot=data_feff_fits_phot[1].data
redshift_phot=Grid_phot.field(0)
weightfn_phot=Grid_phot.field(1)
log10energy_phot=Grid_phot.field(2)
feffbest_phot_o=Grid_phot.field(3)

log10energy_array_phot=zeros(40)
feffbest_array_phot=zeros(40)

for i in range (1,40):
    log10energy_array_phot[i]=log10energy_phot[0,i]

for i in range (1,40):
    feffbest_array_phot[i]=feffbest_phot_o[0,i]

feffbest_phot=interp1d(log10energy_array_phot,feffbest_array_phot)


##     feff for e^+e^-   ##

Species_2="elec"
data_feff_fits_elec=fits.open('feff_summary_'+Species_2+'_'+AssumedIonizationHistory+'.fits')
Grid_elec=data_feff_fits_elec[1].data
redshift_elec=Grid_elec.field(0)
weightfn_elec=Grid_elec.field(1)
log10energy_elec=Grid_elec.field(2)
feffbest_elec_o=Grid_elec.field(3)

log10energy_array_elec=zeros(40)
feffbest_array_elec=zeros(40)

for i in range (1,40):
    log10energy_array_elec[i]=log10energy_elec[0,i]

for i in range (1,40):
    feffbest_array_elec[i]=feffbest_elec_o[0,i]

feffbest_elec=interp1d(log10energy_array_elec,feffbest_array_elec)




#################################################################
##                     Compute the feff                        ##
#################################################################

def feff(mass):
    integrand =lambda logeng: (np.exp(logeng)**2)*2.*Energy_Spectrum_electrons(np.exp(logeng),mass)*feffbest_elec(np.log10(np.exp(logeng)-511000.)) + (np.exp(logeng)**2)*Energy_Spectrum_photons(np.exp(logeng),mass)*feffbest_phot(np.log10(np.exp(logeng)-511000.))
    return (1./(2.*mass))*integrate.quad(integrand,logeng_min(mass),log(mass))[0]

mass=linspace(140.e6,1000.e6,50)
feff_arr=zeros(50)
for i in range(50):
    feff_arr[i]=feff(mass[i])





