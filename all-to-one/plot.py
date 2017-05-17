import json
import numpy as np
import pylab as plt
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

def convert_theta(N, theta):
    W = theta[:N]
    K = np.array(theta[N:]).reshape((N,-1))
    K = array_without_diag(K)
    return W, K

########################################
model_config = 'model.config'
model_file_name = 'model.pkl'

########################################
# Reading model
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

y0 = np.array(params["y0"], dtype=np.float64)
theta = np.array(params["theta"], dtype=np.float64)

########################################
## Reconstruct signal
Y0_rec = y0
W_rec, K_rec = convert_theta(len(y0), theta)

Y0_rec, W_rec, K_rec = Y0_rec[:N], W_rec[:N], K_rec[:N,:N]
y_rec = compute_kuramoto(ts, Y0_rec, W_rec, K_rec)

########################################
## Plot for each oscillator 
for osc in range(N):
    plt.subplot(N, 2, 2*osc+1)
    plt.plot(ts, y_in[osc], 'g')
    plt.plot(ts, y_rec[osc], 'r')
    plt.title("Osc "+str(osc+1))

    plt.subplot(N, 2, 2*osc+2)
    plt.plot(ts, y_in[osc]-y_rec[osc], 'g')
    plt.title("Pointwise diff")

plt.tight_layout()
plt.show()
