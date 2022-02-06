import numpy as np
from scipy.special import legendre
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt


''' Computes the cross section for sound scattering
    by a sphere. The propagation medium is assumed 
    infinite with uniform sound speed.
    
    The wave equation is,
    
        $\nabla^2 \psi + k^2 \psi = 0
'''


# ------------- SETTINGS ----------------- #
c_air = 331.     # speed of sound in air m/s
c_water = 1490.  # speed of sound in water m/s
R = 0.02         # radius of air sphere
l_max = 5       # max orbital angular momentum to consider
freq_min = 100   # Hz
freq_max = 20000 # Hz
plot_only_sum = False
thetas = [180, 0] # detection angles
# ---------------------------------------- #


f_res = c_air / (2 * R)

print('Naive estimates:')
print(f' - Resonance frequency (c/2R):  {f_res:.0f} Hz')
print(f' - Cross section (pi*R^2):      {np.pi*R**2*1e4:.0f} cm^2')


def phase_shift(k, K, R, l):
    """ computes the phase shift using Equation 2.18 in the scattering primer
        k: wave number in outer region, m^-1
        K: wave number in inner region, m^-1
        R: radius of inner region
        l: angular momenutm        
    """
    a = k * spherical_jn(l, K*R) * spherical_jn(l+1, k*R)
    b = K * spherical_jn(l, k*R) * spherical_jn(l+1, K*R)
    c = k * spherical_jn(l, K*R) * spherical_yn(l+1, k*R)
    d = K * spherical_yn(l, k*R) * spherical_jn(l+1, K*R)
    tan_delta_l = (a - b) / (c - d)
    return np.arctan(tan_delta_l)

def cross_section(d, k, l):
    """ angle-integrated cross section in m^2 
        d: phase shift in radians
        k: wave number, m^-1
        l: angular momenutm
    """
    return 4. * np.pi / k**2 * (2*l + 1) * np.sin(d)**2
    
def amplitude(d, k, l, theta):
    """ returns the scattering amplitude at the given angle, f(theta),
        as computed using Messiah X.31 
        d: phase shift
        k: wave number, m^-1
        l: angular momentum
        theta: angle in radians
    """
    Pl = legendre(l)
    f_l = 1. / k * (2*l + 1) * np.exp(1j * d) * np.sin(d) * Pl(np.cos(theta))
    return f_l
    

fig, ax = plt.subplots(2, 1, sharex=True)

freq = np.logspace(np.log10(freq_min), np.log10(freq_max), num=1000)

k = 2*np.pi * freq / c_water
K = 2*np.pi * freq / c_air

# angle-integrated cross section
sigma_tot = 0
for l in range(l_max):
    delta = phase_shift(k=k, K=K, R=R, l=l)
    sigma = cross_section(d=delta, k=k, l=l)
    sigma_tot += sigma
    if not plot_only_sum: ax[0].plot(freq, sigma * 1e4, label=f'l={l}')

ax[0].plot(freq, sigma_tot * 1e4, label=f'sum')

ax[0].legend()

ax[0].set_xlabel('frequency (Hz)')
ax[0].set_ylabel('cross section (cm^2)')


# cross section at specific angles
for theta in thetas:
    f_tot = 0
    for l in range(l_max):
        delta = phase_shift(k=k, K=K, R=R, l=l)
        f_l = amplitude(d=delta, k=k, l=l, theta=theta*np.pi/180.)
        f_tot += f_l

    sigma = np.abs(f_tot)**2

    ax[1].plot(freq, sigma, label=f'{theta:.0f} deg')

ax[1].legend()
ax[1].set_xlabel('frequency (Hz)')
ax[1].set_ylabel(f'cross section (cm^2 / sr)')


plt.show()

