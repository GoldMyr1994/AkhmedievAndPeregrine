import numpy as np
import scipy.fftpack as fftpack
import pandas as pd

################
# AKHMEDIEV
################


def akhmediev_Omega_2_params(Omega):
    a = (1./2)*(1-(Omega**2)/4)
    b = np.sqrt(8*a*(1-2*a))
    return Omega, a, b


def akhmediev_a_2_params(a):
    Omega = 2*np.sqrt(-(2.*a-1.))
    return akhmediev_Omega_2_params(Omega)


def akhmediev(z, t, Omega, show_params=False):

    if type(z) is int or type(z) is float:
        z = np.array([z])

    if type(z) is list:
        z = np.array(z)

    Z, T = np.meshgrid(z, t)

    Omega, a, b = akhmediev_Omega_2_params(Omega)
    if show_params:
        print("akhm: (Omega, a, b)=(%.2f,%.2f,%.2f)" % (Omega, a, b))
    N = (1-4*a)*np.cosh(b*Z)+np.sqrt(2*a)*np.cos(Omega*T)+1j*b*np.sinh(b*Z)
    D = np.sqrt(2*a)*np.cos(Omega*T)-np.cosh(b*Z)
    AKHM = ((N/D)*np.exp(1j*Z)).T
    if AKHM.shape[0] is 1:
        AKHM = AKHM[0]
    return AKHM


################
# PEREGRINE
################

def peregrine(z, t):

    if type(z) is int or type(z) is float:
        z = np.array([z])

    if type(z) is list:
        z = np.array(z)

    Z, T = np.meshgrid(z, t)
    PER = ((1-4*(1+2*1j*Z)/(1+4*T**2+4*Z**2))*np.exp(1j*Z)).T
    if PER.shape[0] is 1:
        PER = PER[0]
    return PER


################
# SIMULATOR
################


def simulate(initial_condition, beta2, gamma, z, t, f, dz, dt):

    # spectra of the initial condition
    initial_condition_spectra = dt*fftpack.fft(initial_condition)

    # list to store the values in time for each z step
    zt_values = [initial_condition]
    # list to store the spectra of values in time for each z step
    zf_values = [initial_condition_spectra]

    # values in time at z=z0
    values = zt_values[0]
    # spectra of values in time at z=z0
    values_spectra = zf_values[0]

    # for each z step
    for _ in z[1:]:

        # LINEAR DISPERSIVE STEP
        _prop = .5*beta2*(2*np.pi*f)**2
        _fact = 1j*_prop*dz
        values_spectra = values_spectra*np.exp(_fact)
        values = (1/dt)*fftpack.ifft(values_spectra)

        # NON LINEAR CHI3 EFFECTS
        values = values*np.exp(1j*gamma*np.abs(values)**2*dz)
        values_spectra = dt*fftpack.fft(values)

        # store results of current step
        zt_values.append(values)
        zf_values.append(values_spectra)

    # dim0, dim1: z, t
    zt_values = np.array(zt_values)
    zf_values = np.array(zf_values)

    return zt_values, zf_values
