import json
import numpy as np
import pylab as plt
import pickle
from scipy.integrate import ode

def kuramoto_ODE(t, y, arg):

    w, k = arg
    yt = y[:,None]
    dy = y-yt
    phase = np.array(w)
   # for m, _k in enumerate(k):
   #     phase += np.sum(_k*np.sin((m+1)*dy),axis=1)
    phase += np.sum(k*np.sin(dy),axis=1)

    return phase

def compute_kuramoto(t, Y0, W, K):
    kODE = ode(kuramoto_ODE)
    kODE.set_integrator("dopri45")
    #kODE.set_integrator("lsoda", nsteps=5000)

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

def convert_theta(theta):
    N = int(np.sqrt(len(theta)+0.25)-0.5)
    W = theta[:N]
    K = np.array(theta[N:]).reshape((N,N))
    K -= np.diag(np.diag(K))
    return W, K

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
W = np.array(config["W"])[:N]
Y = np.array(config["Y"])[:N]
K = np.array(config["K"])[:N,:N]
tMin = config["T"]["tMin"]
tMax = config["T"]["tMax"]
tN = config["T"]["tN"]

## Convert items
K = np.array(K)
W = np.array(W, dtype=np.float64)

ts = np.linspace(tMin, tMax, tN)
y_in = compute_kuramoto(ts, Y, W, K)

########################################
## Load reconstructed signals
with open("op.pkl","rb") as f:
    params = pickle.load(f)

Y_rec = params["y0"]
W_rec, K_rec = convert_theta(params["theta"])
Y_rec, W_rec, K_rec = Y_rec[:N], W_rec[:N], K_rec[:N,:N]
#W_rec = W_rec.astype(np.float64)
W_rec = [float(w) for w in W_rec]
print("Y_rec: ", Y_rec)
print("W_rec: ", W_rec)
print("K_rec: ", K_rec)
y_rec = compute_kuramoto(ts, Y_rec, W_rec, K_rec)
#y_rec = compute_kuramoto(ts, Y_rec, W_rec, K)

########################################
## Plot for each oscillator 
for osc in range(N):
    plt.subplot(N, 2, 2*osc+1)
    plt.plot(ts, y_in[osc], 'g')
    plt.plot(ts, y_rec[osc], 'r')

    plt.subplot(N, 2, 2*osc+2)
    plt.plot(ts, y_in[osc]-y_rec[osc], 'g')
plt.show()
