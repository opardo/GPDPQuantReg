# GPDPQuantReg

## Author

Carlos Omar Pardo (omarpardog@gmail.com)

## Description

R Package. Bayesian and nonparametric quantile regression, using Gaussian Processes to model the trend, and Dirichlet Processes for the error. An MCMC algorithm works behind to fit the model.

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
predictive_data <- data_frame(x = seq(-15, 15, 0.25))
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
