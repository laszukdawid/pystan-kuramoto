functions {

    int to_int(real v){
        int n=0;

        int s = v>=0 ? 1 : -1;
        real abs_v = abs(v);

        while(n<abs_v) n = n + s;

        return n;
    }

    real[] osc(real t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[] x_i) {

        int s = size(theta);
        int N = to_int(sqrt(s+0.25)-0.5);

        real dydt[N];
        real omega[N];
        real K[N,N];

        for (n in 1:N){
            omega[n] = theta[n];

            // Extract K
            for (m in 1:N){
                int idx = (n-1)*N + m;
                K[m][n] = theta[N + idx];
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
    real y_in[T,N];
    real t0;
    real ts[T];
    real<lower=0.001> sig_err;
}
transformed data {
    real x_r[0];
    int x_i[0];
}
parameters {
    real y0[N];
    real theta[N*(1+N)];
}
transformed parameters {
}
model {

    real y_hat[T,N];

    y_hat = integrate_ode_rk45(osc, y0, t0, ts, theta, x_r, x_i);
    for (t in 1:T){
        for (n in 1:N){
            y_in[t,n] ~ normal(y_hat[t,n], sig_err);
        }
    }
}
generated quantities {
}
