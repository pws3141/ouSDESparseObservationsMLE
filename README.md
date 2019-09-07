# Maximum Likelihood Estimator for Sparse Observations

Assuming time series Y_t is of the form:
> dX_t = - \gamma X_t dt + \sigma dB_t;
> Y_t = X_t + m(t),

where m(t) is estimated using quadratic splines.

Maximum likelihood estimator is obtained by optimising the log-likelihood found from
considering the transition density between each point, which is available in
closed-form for the OU process X_t.

## Spline estimation

To obtain m(t), use [Spline Estimation](200~https://github.com/pws3141/splineEstimation)
function.

