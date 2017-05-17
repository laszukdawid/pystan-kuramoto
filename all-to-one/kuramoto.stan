functions {

    real[] osc(real t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[] x_i) {

        int N = x_i[1];

        real dydt[N];
        real omega[N];
        real K[N,N];

        int idx = 0;
        for (n in 1:N){
            omega[n] = theta[n];

            // Extract K
            for (m in 1:N){
                if (m==n) {
                    K[m][n] = 0;
                } else {
                    idx = idx + 1;
                    K[m][n] = theta[N + idx];
                }
            }
        }

        for (n in 1:N){
            dydt[n] = omega[n];

            for (m in 1:N){
                dydt[n] = dydt[n] + K[m][n]*sin(y[m]-y[n]);
            }
        }

        return dydt;
    }
}
data {
    int<lower=1> T;
    int N;
    real y_in[T];
    real t0;
    real ts[T];
    real<lower=0.001> sig_err;
}
transformed data {
    real x_r[0];
    int x_i[1];

    x_i[1] = N;
}
parameters {
    real y0[N];
    real theta[N*N];
}
transformed parameters {
}
model {

    real y_hat[T,N];
    real y_tmp;

    y_hat = integrate_ode_rk45(osc, y0, t0, ts, theta, x_r, x_i);
    for (t in 1:T){
        y_tmp = 0;
        for (n in 1:N){
            y_tmp = y_tmp + y_hat[t,n];
        }
        y_in[t] ~ normal(y_tmp, sig_err);
    }
}
generated quantities {
}
