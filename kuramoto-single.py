#/usr/bin/python
from __future__ import division, print_function
import json
import numpy as np
import pystan
import pickle
from scipy.integrate import ode

def array_without_diag(arr):
    tmp_arr = np.eye(arr.shape[0])
    tmp_arr[tmp_arr!=1] = arr.flatten()
    np.fill_diagonal(tmp_arr, 0)
    return tmp_arr

def kuramoto_ODE(t, y, arg):

    w, k = arg
    yt = y[:,None]
    dy = y-yt
    phase = np.array(w)
    phase += np.sum(k*np.sin(dy),axis=1)

    return phase

def compute_kuramoto(t, Y0, W, K):
    # Assuming two options: square or lacking diag elems
    print(K)
    if len(K.shape)>1 and K.shape[0]!=K.shape[1]:
        K = array_without_diag(K)

    kODE = ode(kuramoto_ODE)
    kODE.set_integrator("dopri5")

    # Set parameters into model
    kODE.set_initial_value(Y0, t[0])
    kODE.set_f_params((W, K))

    phase = np.empty((len(W), len(t)))

    # Run ODE integrator
    for idx, _t in enumerate(t[1:]):
        phase[:,idx] = kODE.y
        kODE.integrate(_t)

    phase[:,-1] = kODE.y
    return phase

########################################
model_stan = 'kuramoto.stan'
model_file_name = 'model.pkl'
sample_file = 'samples'
model_config = 'model.config'

########################################
# Reading model
with open(model_stan,'r') as f:
    ocode = f.read()

with open(model_config, 'r') as f:
    config = json.loads(f.read())
    print(config)

########################################
########## EXTRACT PARAMS
print("Defining parameters")
N = config["N"]
W = np.array(config["W"], dtype=np.float64)[:N]
Y = np.array(config["Y"], dtype=np.float64)[:N]
K = np.array(config["K"], dtype=np.float64)[:N,:N]
tMin = config["T"]["tMin"]
tMax = config["T"]["tMax"]
tN = config["T"]["tN"]

## Convert items
np.fill_diagonal(K, np.inf)
K_flat = K.flatten()
K_flat = K_flat[K_flat!=np.inf]
np.fill_diagonal(K, 0)

########################################
print("Starting model...", end='')
try:
    model = pickle.load(open(model_file_name, 'rb'))
    print(" managed to load it")
except Exception as e:
    print(" needs compilation")
    print(e)
    model = pystan.StanModel(model_code=ocode, verbose=True)
    model.show()

    print("Saving results to a file: " + model_file_name)
    with open(model_file_name, 'wb') as f:
        pickle.dump(model, f)

########################################
## Compute

ts = np.linspace(tMin, tMax, tN)
y_in = compute_kuramoto(ts, Y, W, K)
y_in = np.sum(y_in, axis=0)

theta = np.append(W, K_flat)
init_theta = np.append(Y,theta)
print("init_theta: ", init_theta)

data = dict(T=tN, N=N, ts=ts, y_in=y_in, t0=tMin-0.0001,  sig_err=.05)
data["init_params"] = init_theta
init = dict(theta=init_theta)
print("Optimising with data")
rand = lambda n: np.random.random(n)

init_param = {"y0":Y+rand(len(Y)), "theta": theta+rand(len(theta))}
#init_param = {"y0":Y, "theta": theta}

op = model.optimizing(data=data, init=init_param, sample_file=sample_file)

print("Presenting results")
print(op)


with open('op.pkl', 'wb') as f:
    pickle.dump(op, f)
