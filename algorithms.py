import numpy as np
import math
import random

from config import ArgProcessor

arguments = ArgProcessor()
input_data = arguments.getInputs()

DEFAULT_AMPLITUDE = input_data['amplitude'] #0.6 ~ mean for best model
WIND_SPEED = input_data['wind_speed'] #1.5 ~ mean for best model
MIN_WAVE_SIZE = input_data['min_wave_size'] #0.006 ~ mean for best model, magnitude -> min => wave => max
SURFACE_DIM = input_data['grid_dim'] #100 ~ mean for best model

CARTESIAN_POINTS_DIMENSION = 200

LX = 20.0
LY = 20.0
THETA = 45
X_0 = 0.1

#L = V^2/g
#each wave component mass/density
m_wave = 0.01
m = m_wave*CARTESIAN_POINTS_DIMENSION
g = 9.82

# h ~ Amplitude of waves ~ <A_i>
h_wave = 10.0
h_deep = 5.0
_lambda = 1.5

c = math.sqrt(g * h_deep)
k = 2.0*math.pi/_lambda

omega_deep = np.sqrt(g*k)
omega_shallow = np.sqrt(g*k*np.tanh(k*h_deep))

#Gaussian random generator.
#The numbers are generated
#using the Box-Muller method.
#Also we can generate mu, sigma by pseudo random values generator
#if want to use classic model.
def Gaussian():

    while True:
        p1 = random.uniform(0.0, 1.0)
        p2 = random.uniform(0.0, 1.0)

        z1 = np.sqrt(-2.0 * np.log(p1)) * np.cos(2.0 * np.pi * p2)
        z2 = np.sqrt(-2.0 * np.log(p1)) * np.sin(2.0 * np.pi * p2)
        result_d = np.round((z1 + z2)/2.0, 4)

        if result_d >= 1.0:
            return result_d

#Direct Fourier Transform
#approximate solution, get new phases - shiftness waves factor
def DFT(x):
    N = x.size
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * np.pi * k * n / N)
    return np.dot(e, x)

#Gerstner OX wave coordinate component
def OXCartesianGerstner(t, A = DEFAULT_AMPLITUDE):
    global DEFAULT_AMPLITUDE
    x = X_0 - k*A*np.sin(k*X_0 - omega_deep*t)
    return x

#Gerstner OY wave coordinate component
def OYCartesianGerstner(t, A = DEFAULT_AMPLITUDE):
    global DEFAULT_AMPLITUDE
    y = A*np.cos(k*X_0 - omega_shallow*t)
    return y

#Gerstner OX wave velocity component
def VXCartesianGerstner(_coordinates_x):
    dx = (np.max(_coordinates_x) - np.min(_coordinates_x))/_coordinates_x.size
    grad_x = np.gradient(_coordinates_x, dx)
    return grad_x

#Gerstner OY wave velocity component
def VYCartesianGerstner(_coordinates_y):
    dy = (np.max(_coordinates_y) - np.min(_coordinates_y))/_coordinates_y.size
    grad_y = np.gradient(_coordinates_y, dy)
    return grad_y

#Product factor
def Grid(a, b):
    return a*b

#Philipps potential - part of statistical distribution of wave/big-wave amplitudes
#Stochastic factor raised by product with gaussian magnitude c <[1, 2]>
def PhillipsPotential(x, y):
    global DEFAULT_AMPLITUDE
    kx = (2.0*np.pi*x)/2.0
    ky = (2.0*np.pi*y)/2.0
    k_sq = kx*kx + ky*ky
    L_sq = np.power(WIND_SPEED*WIND_SPEED/g, 2.0)

    if k_sq.any() == 0.0:
        return 0.0

    wind_alignment = 2.0

    #if we want use to more precise amplitude distribution
    #but this value, for linear modelling exist in a small range of values
    #and correlating as: a -> max, wave size -> min, otherwise it starts nonlinear process
    # A = np.sqrt(g * np.sqrt(np.power((2.0 * np.pi * x) / LX, 2.0) + np.power((2.0 * np.pi * y) / LY, 2.0)) *
    #             (1.0 + (np.power((2.0 * np.pi * x) / LX, 2) + np.power((2.0 * np.pi * y) / LY, 2.0)) * np.power(L, 2.0)))

    A = DEFAULT_AMPLITUDE

    potential =  (A*np.exp((-1.0/(k_sq*L_sq))) * np.exp(-k_sq*np.power(MIN_WAVE_SIZE, 2.0)) * np.power((kx*kx)/k_sq, wind_alignment))/k_sq*k_sq

    # print('[kx = {}\tky = {}\tk_sq = {}\tL_sq = {}\tA = {}]'.format(np.round(kx, 2), np.round(ky, 2), np.round(k_sq, 2), np.round(L_sq, 2), np.round(A, 2)))
    # print('A = {} B = {} C = {} D = {}'.format(A*np.exp((-1.0/(k_sq*L_sq))), np.exp(-k_sq*np.power(MIN_WAVE_SIZE, 2.0)), np.power((kx*kx)/k_sq, wind_alignment), k_sq*k_sq))
    # print('P(x, y) = {}'.format(np.round(potential, 2)))
    return potential

#bspline series
#ortogonal projection - u (x) or v (y) waves
def b_n_series(dx, dy, t = 1.0):

    omega_product = 1.0
    b_u = 0.0
    bn = 0.0

    _du = 0.5
    for i in np.arange(0.1, 1.0, _du):

        u_i = dx
        u_j = dy

        cron_u = 1.0

        if np.sqrt(u_i**2.0 + u_j**2.0) >= u_i:
            cron_u = np.power((np.sqrt(u_i**2.0 + u_j**2.0) - u_i), _du)

        if np.sqrt(u_i**2.0 + u_j**2.0) < u_i:
            cron_u = 0.0

        omega_product *= (1.0/(u_j - u_i))*cron_u*np.sign(u_j)

        b_u += (i+1.0/i)*omega_product

        bn += omega_product*b_u*t

        #print('ui = {}'.format(u_i))
        #print('uj = {}'.format(u_j))
        #print('bn = {}'.format(bn))

    return bn

#translate generated points
#we can put there maximum of algorithm logic and modify x, y, step, while processing
def CartesianDistibutionpoints(_step, _dim = CARTESIAN_POINTS_DIMENSION):
    f_ox = list()
    f_oy = list()
    _du = _step

    for i in np.arange(_du, _du*_dim + _du, _du):
        f_ox.append(OXCartesianGerstner(i))
        f_oy.append(OYCartesianGerstner(i))

    f_ox = np.array(f_ox)
    f_oy = np.array(f_oy)

    return [f_ox, f_oy]

#translate x cartesian coordinate to u cylindrical coordinate
def CartesianOXtoCylindricalU(x, y):
    _r = np.sqrt(x**2.0 + y**2.0)
    _phi = np.arctan2(y, x)
    return _r*np.cos(_phi)

#translate y cartesian coordinate to v cylindrical coordinate
def CartesianOYtoCylindricalV(x, y):
    _r = np.sqrt(x**2.0 + y**2.0)
    _phi = np.arctan2(y, x)
    return _r*np.sin(_phi)

#OZ projection in cylindrical basis equal to cartesian
