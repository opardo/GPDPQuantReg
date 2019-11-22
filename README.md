# GPDPQuantReg

## Author

Carlos Omar Pardo Gomez (cop2108@columbia.edu)

## Overview

R Package. Bayesian and nonparametric quantile regression, using Gaussian Processes to model the trend, and Dirichlet Processes for the error. An MCMC algorithm works behind to fit the model.

## Model

Be one dependent random variable _y_, and a vector of independent variables _x_, so _P(y|x)_ makes sense. Be _q<sub>p</sub>_ the function which returns the _p_-quantile for a random variable. The model assumes that given any probability _p_ (between 0 and 1) for which you want to estimate the _y|x_'s _p_-quantile, an observation _y_ comes from the sum

<a href="https://www.codecogs.com/eqnedit.php?latex=y&space;=&space;f_p(x)&space;&plus;&space;\varepsilon_p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y&space;=&space;f_p(x)&space;&plus;&space;\varepsilon_p" title="y = f_p(x) + \varepsilon_p" /></a>,

where _f<sub>p</sub>_ is the function which links _x_ to _y_, and epsilon is a dispersion random variable, such that its _p_-quantile is 0. Then, because the quantile of a sum is the sum of the quantiles, we get that  _q<sub>p</sub>(y|x)_ = _f<sub>p</sub>(x)_. This package is going to focus on estimating _f<sub>p</sub>_.

We are going to address that problem by using the __GPDP model__, which is described below. (For a complete understanding of it, it is recommended to have some previous knowledge in _Gaussian Processes_ ([Bagnell's lecture](http://www.cs.cmu.edu/~16831-f14/notes/F09/lec21/16831_lecture21.sross.pdf)) and _Dirichlet Processes_ ([Teh 2010](https://www.stats.ox.ac.uk/~teh/research/npbayes/Teh2010a.pdf)), particularly, the stick-breaking representation). 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\dpi{150}&space;$$\begin{aligned}&space;y_i|&space;f_p(x_i),&space;z_i,&space;\sigma_k^*&space;&\sim&space;AL_p({\varepsilon_p}_i&space;=&space;y_i&space;-&space;f_p(x_i)&space;|&space;\sigma_{z_i}),&space;\\&space;f_p|m,&space;k,&space;\lambda&space;&\sim&space;\mathcal{GP}(m,k(\lambda)|\lambda),&space;\\&space;\lambda&space;&\sim&space;GI(c_\lambda,d_\lambda),&space;\\&space;z_i&space;|&space;\pi&space;&\sim&space;Mult_\infty(\pi),&space;\\&space;\pi&space;|&space;\alpha&space;&\sim&space;GEM(\alpha),&space;\\&space;\sigma_k^*&space;|&space;c_{DP},&space;d_{DP}&space;&\sim&space;GI(\sigma_k|c_{DP},&space;d_{DP}),\\&space;k(x_i,&space;x_j&space;|&space;\lambda)&space;&=&space;\lambda&space;\text{&space;}&space;exp\{-\norm{x_i&space;-&space;x_j}_2\}.&space;\end{aligned}$$" target="_blank"><img src="https://latex.codecogs.com/png.latex?\inline&space;\dpi{150}&space;$$\begin{aligned}&space;y_i|&space;f_p(x_i),&space;z_i,&space;\sigma_k^*&space;&\sim&space;AL_p({\varepsilon_p}_i&space;=&space;y_i&space;-&space;f_p(x_i)&space;|&space;\sigma_{z_i}),&space;\\&space;f_p|m,&space;k,&space;\lambda&space;&\sim&space;\mathcal{GP}(m,k(\lambda)|\lambda),&space;\\&space;\lambda&space;&\sim&space;GI(c_\lambda,d_\lambda),&space;\\&space;z_i&space;|&space;\pi&space;&\sim&space;Mult_\infty(\pi),&space;\\&space;\pi&space;|&space;\alpha&space;&\sim&space;GEM(\alpha),&space;\\&space;\sigma_k^*&space;|&space;c_{DP},&space;d_{DP}&space;&\sim&space;GI(\sigma_k|c_{DP},&space;d_{DP}),\\&space;k(x_i,&space;x_j&space;|&space;\lambda)&space;&=&space;\lambda&space;\text{&space;}&space;exp\{-||x_i&space;-&space;x_j||_2\}.&space;\end{aligned}$$" title="$$\begin{aligned} y_i| f_p(x_i), z_i, \sigma_k^* &\sim AL_p({\varepsilon_p}_i = y_i - f_p(x_i) | \sigma_{z_i}), \\ f_p|m, k, \lambda &\sim \mathcal{GP}(m,k(\lambda)|\lambda), \\ \lambda &\sim GI(c_\lambda,d_\lambda), \\ z_i | \pi &\sim Mult_\infty(\pi), \\ \pi | \alpha &\sim GEM(\alpha), \\ \sigma_k^* | c_{DP}, d_{DP} &\sim GI(\sigma_k|c_{DP}, d_{DP}),\\ k(x_i, x_j | \lambda) &= \lambda \text{ } exp\{-||x_i - x_j||_2\}. \end{aligned}$$" /></a>

Where:
- _p_ is the probability for which you want to estimate the quantile.
- _AL<sub>p</sub>_ is the Assymetric Laplace Distribution (ALD), with density function given by <a href="https://www.codecogs.com/eqnedit.php?latex=w_p^{AL}(u|\sigma)&space;=&space;\frac{p(1-p)}{\sigma}&space;exp\left[&space;-\rho_p&space;\left(&space;\frac{u}{\sigma}&space;\right)&space;\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?w_p^{AL}(u|\sigma)&space;=&space;\frac{p(1-p)}{\sigma}&space;exp\left[&space;-\rho_p&space;\left(&space;\frac{u}{\sigma}&space;\right)&space;\right]" title="w_p^{AL}(u|\sigma) = \frac{p(1-p)}{\sigma} exp\left[ -\rho_p \left( \frac{u}{\sigma} \right) \right]" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=\rho_p(u)&space;=&space;u&space;\times&space;[pI_{(u>0)}&space;-&space;(1-p)&space;I_{(u<0)})]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho_p(u)&space;=&space;u&space;\times&space;[pI_{(u>0)}&space;-&space;(1-p)&space;I_{(u<0)})]" title="\rho_p(u) = u \times [pI_{(u>0)} - (1-p) I_{(u<0)})]" /></a>.
- _GP_ is a Gaussian process, with mean (_m_) and covariance (_k_) functions. 
- _GI_ is the Inverse-Gamma distribution.
- _Mult<sub>inf</sub>_ is the Multinomial distribution, when the number of categories tend to infinity.
- GEM (for Griffiths, Engen and McCloskey) is a distribution used in Dirichlet Processes' literature, as described in Teh (2010).

## Algorithm

An MCMC algorithm is used to find the _f<sub>p</sub>_'s posterior distribution via simulations, particularly, a _Gibbs sampler_ has been developed.

Since the theoretical model is a nonparametric one, it contemplates infinite parameters, particularly for the Dirichlet process. However, it's clear we cannot estimate and allocate such a number of values, so this packages uses the _slice sampling_ algorithm proposed by [Kalli et al.](http://users.wpi.edu/~balnan/Kalli-Griffin-Walker-2011.pdf) to truncate the number of them in a dynamic way. The results approximately converge to the expected ones, if we could do it in the theoretical way.

## Example

First, you must install the package directly from Github. (I hope eventually you will do it from CRAN.)

```{r}
library(devtools)
install_github("opardo/GPDPQuantReg")
library(GPDPQuantReg)
```

Then, you can create some artificial data to start familiarizing with the package's dynamic. In this case, 
I'm using a complex trend function and a non-normal error, sampling ONLY 20 POINTS.

```{r}
set.seed(201707)
f_x <- function(x) return(0.5 * x * cos(x) - exp(0.1 * x))
error <- function(m) rgamma(m, 2, 1)
m <- 20
x <- sort(sample(seq(-15, 15, 0.005), m))
sample_data <- data.frame(x = x, y = f_x(x) + error(m))
```

![Image of sample_data](https://github.com/opardo/GPDPQuantReg/blob/master/images/sample_data.png)

Now, it's time to fit the model with a MCMC algorithm for a specific _p_ probability.

```{r}
GPDP_MCMC <- GPDPQuantReg(y ~ x, sample_data, p = 0.250)
```

Since it is a nonparametric model, it is focused on prediction with a credible interval.

```{r}
predictive_data <- data.frame(x = seq(-15, 15, 0.25))
credibility <- 0.90
prediction <- predict(GPDP_MCMC, predictive_data, credibility)
```
And for a complex trend function, non-normal error and only 20 sampled points... we get AWESOME RESULTS!

![Image of prediction](https://github.com/opardo/GPDPQuantReg/blob/master/images/prediction.png)

Some diagnostics (ergodicity, autocorrelation, crosscorrelation and traces) for the Markov Chains are available too.

```{r}
diagnose(GPDP_MCMC)
```

![Image of ergodicity](https://github.com/opardo/GPDPQuantReg/blob/master/images/ergodicity.png)
![Image of autocorrelation](https://github.com/opardo/GPDPQuantReg/blob/master/images/autocorrelation.png)
![Image of crosscorrelation](https://github.com/opardo/GPDPQuantReg/blob/master/images/crosscorrelation.png)
![Image of trace](https://github.com/opardo/GPDPQuantReg/blob/master/images/trace.png)

## Test

TODO
