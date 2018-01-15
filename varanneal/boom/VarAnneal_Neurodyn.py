
import numpy as np
from scipy import interpolate
from varanneal import va_ode
import os, sys, time
import pickle

import matplotlib
matplotlib.use('nbagg')
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

from matplotlib import gridspec

# For 3D plots
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


from varanneal import va_ode

# For running this on the cluster
initID = int(sys.argv[1])
adolcID = int(sys.argv[2])

# ## Define the ODE system

# Load in the input current


Ipath = os.getcwd() + "/IforRealNeuron.csv"
Idat = np.genfromtxt(Ipath, delimiter=',')

# Load in the parameters
sf = np.load("model_parameters/scale_factors.npy")
with open("model_parameters/model_params.txt", "rb") as fp:
    model_params = pickle.load(fp)

# Unpack scale factors
V_scale = sf[0]
C_scale = sf[1]
I_scale = sf[2]
I_inj_scale = sf[3]
R_scale = sf[4]
# Unpack parameters
mp = model_params
C_m = 4e-12 * C_scale
g = mp[0]
E_rev = mp[1]
vBias = mp[2]
am = mp[3]
bm = mp[4]
ah = mp[5]
bh = mp[6]
an = mp[7]
bn = mp[8]
# should add this
I_sampling_rate = 5e3


# Define alpha and beta representations
def sigma(vBiask, V, sign = 1):
    mu = 0.7
    Ut = 26e-3 * V_scale
    return 1 / (1 + np.exp(sign * mu * (vBiask - V) / Ut))


def alpha_spline(V, x, vBias=vBias, am=am, ah=ah, an=an):
    """
    Used to compute the conductance at each time point. 
    The optional arguments default to the parameters defined at the start of the document
    """
    alpha = 0
    for k in np.arange(7):
        if x == "m":
            alpha += am[k] * sigma(vBias[k], V, 1)
        if x == "h":
            alpha += ah[k] * sigma(vBias[k], V, -1)
        if x == "n":
            alpha += an[k] * sigma(vBias[k], V, 1)
    return alpha


def beta_spline(V, x, vBias=vBias, bm=bm, bh=bh, bn=bn):
    beta = 0
    for k in np.arange(7):
        if x == "m":
            beta += bm[k] * sigma(vBias[k], V, -1)
        if x == "h":
            beta += bh[k] * sigma(vBias[k], V, 1)
        if x == "n":
            beta += bn[k] * sigma(vBias[k], V, -1)
    return beta


# Prepare injected current


#Used to interpolate time points that are undefined in Idat
fIdat = interpolate.interp1d(np.arange(0,len(Idat)), Idat) 

# should save the sampling rate and load it in
def I_inj(t, scale=I_inj_scale):
    if t[0] * I_sampling_rate <= len(Idat):
        return fIdat(t * I_sampling_rate) * scale
    else:
        return 0

# Define model


def neuron(t, y, k):
    #v, m, h, n = y
    v = y[:,0]
    m = y[:,1]
    h = y[:,2]
    n = y[:,3]
    # g = (2.62e-8, 5.25e-9, 4.9e-10)
    # E_rev = (1.109, 0.923, 0.9304)
    g = (k[0], k[1], k[2])
    E_rev = (k[3], k[4], k[5])
    
    I_na = g[0] * m**3 * h * (v - E_rev[0])
    I_k = g[1] * n**4 * (v - E_rev[1])
    I_l = g[2] * (v - E_rev[2])

    dvdt = (I_inj(t) - I_na - I_l - I_k) / C_m
    dmdt = alpha_spline(v, "m") * (1 - m) - beta_spline(v, "m") * m
    dhdt = alpha_spline(v, "h") * (1 - h) - beta_spline(v, "h") * h
    dndt = alpha_spline(v, "n") * (1 - n) - beta_spline(v, "n") * n
    dydt = np.transpose(np.array([dvdt, dmdt, dhdt, dndt]))
    
    return dydt




def nakl(t, y, P):
    """
    Neuron Model
    
    Paul's example has P as Pstim, which includes the stimulating current. Not sure if that is more efficient
    """
    v, m, h, n = (y[:,0], y[:,1], y[:,2], y[:,3])
    # Load parameters
    g = (P[0], P[1], P[2])
    E_rev = (P[3], P[4], P[5])
    vBias = (P[6:13])
    
    am = (P[13:20])
    bm = (P[20:27])
    ah = (P[27:34])
    bh = (P[34:41])
    an = (P[41:48])
    bn = (P[48:55])
    
    I_na = g[0] * m**3 * h * (v - E_rev[0])
    I_k = g[1] * n**4 * (v - E_rev[1])
    I_l = g[2] * (v - E_rev[2])

    dydt = np.zeros_like(y)

    dydt[:,0] = (I_inj(t) - I_na - I_l - I_k) / C_m
    dydt[:,1] = alpha_spline(v, "m", vBias) * (1 - m) - beta_spline(v, "m", vBias) * m
    dydt[:,2] = alpha_spline(v, "h", vBias) * (1 - h) - beta_spline(v, "h", vBias) * h
    dydt[:,3] = alpha_spline(v, "n", vBias) * (1 - n) - beta_spline(v, "n", vBias) * n
    
    return dydt


# #### Action/annealing (hyper)parameters


# Model system dimension
D = 4

# Measured variable indices
# (-, t) (0, v) (1, m) (2, h) (3, n) (4, I)
Lidx = [0, 1, 2, 3]

# RM, RF0
RM = 1.0 / (0.5**2)
RF0 = 4.0e-6

# alpha, and beta ladder
alpha = 1.3
beta_array = np.linspace(0, 200, 201)
#beta_array = np.linspace(0, 10, 11)

g0 = RF0/RM
gammas_all = g0 * alpha**beta_array


# #### Load observed data



data = np.load(os.getcwd() + "/ode_data.npy")
times_data = data[:, 0]
dt_data = times_data[1] - times_data[0]
N_data = len(times_data)

#extracting observed data here
data = data[:, 1:]
data = data[:, Lidx]


# Set $\Delta t_f$ based on $\Delta t$.


# model state discretization
freq_mod = 1.0  # how often to put down a state variable
dt_model = dt_data / freq_mod
if freq_mod == 1.0:
    N_model = N_data
else:
    N_model = int(N_data * freq_mod) - 1


# #### Initial path/parameter guesses
# Later in the notebook, we'll have the option of setting the initial guesses for the observed variables equal to the observations themselves.


# State variables
# This should be an array with N_f elements, where element n_f is a D-dimensional 
# vector. In other words, this is an array of state vectors at all "model times".
X0 = (20.0*np.random.rand(N_model * D) - 10.0).reshape((N_model, D))

#n_params = 6
n_params = 55
# Parameters
Pidx = np.arange(0,n_params)


# Setting ranges for initial guesses
# r controls the size of the search space of the parameters
r = 0.5
Pg = []

# In case we want to estimate a subset of the full parameter set
param_count = 0

for grp, params in enumerate(model_params):
    for p in params:
        param_count += 1
        Pg.append([p * (1 - r), p * (1 + r)])
    if param_count == n_params:
        break


# Compiling initial guesses
Pinit = np.zeros(len(Pidx))

# seed it!
np.random.seed(3918420 + initID)

# Initial guesses
for i, b in enumerate(Pg):
    r = b[1] - b[0]
    Pinit[i] = r*np.random.rand() + b[0]

# cast it as a numpy array
Pinit = np.array(Pinit)


# #### Use VA to estimate states and parameters
# First we need to initialize an Annealer object, which stores information about the model, data, annealing hyperparameters, and the action.  It also executes the VA algorithm, then is used to save the state and parameter estimates to file.

# Initialize Annealer
anneal1 = va_ode.Annealer()

# Set the model
if n_params == 6:
    anneal1.set_model(neuron, D)
else:
    anneal1.set_model(nakl, D)

# Load the data into the Annealer object
anneal1.set_data(data, t=times_data)


# Run VA


# First set some options for the optimization.
# Bounds [v], [m], [h], [n]
bounds = [[0 * V_scale, 1.5 * V_scale], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
for i in xrange(len(Pidx)):
    bounds.append(Pg[i])
    
    
# The full list of options can be found in the scipy.optimization package documentation.
BFGS_options = {'gtol':1.0e-8, 'ftol':1.0e-8, 'maxfun':1000000, 'maxiter':1000000}

tstart = time.time()  # time how long VA takes

# Annealer.anneal() executes VA for all beta values (defined above)
# Note the init_to_data option: this initializes the measured variables to the data.
anneal1.anneal(X0, Pinit, alpha, beta_array, RM, RF0, Lidx, Pidx, dt_model=dt_model,
               init_to_data=True, disc='SimpsonHermite', method='L-BFGS-B',
               opt_args=BFGS_options, adolcID=adolcID)

print("\nADOL-C annealing completed in %f s."%(time.time() - tstart))


# # Save action, constituent errors, and state/parameter estimates to file.


anneal1.save_paths("results/paths/paths_%d.npy" % (initID)) #state paths
anneal1.save_params("results/params/params_%d.npy" % (initID))
anneal1.save_action_errors("results/boom/action_errors/action_errors%d.npy" % (initID))#saves action and constituent errors


