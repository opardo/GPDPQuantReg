# GPDPQuantReg

## Author

Omar Pardo (omarpardog@gmail.com)

## Description

Bayesian and nonparametric quantile regression, using Gaussian Processes to model the base function and Dirichlet Processes
for the error, and a MCMC algorithm to fit the model.

## Example

First, you must install the package directly from Github.

```{r}
library(devtools)
install_github("opardo/GPDPQuantReg")
library(GPDPQuantReg)
```

Then, you can create some artificial data to start familiarizing with the package's dynamic. In this case, 
I'm using a complex base function and a non-normal error.

```{r}
set.seed(201707)
f_x <- function(x) return(0.5 * x * cos(x) - exp(0.1 * x))
error <- function(m) rgamma(m, 1, 1)
m <- 20
x <- sort(sample(seq(-15, 15, 0.005), m))
sample_data <- data.frame(x = x, y = f_x(x) + error(m))
```

Now, it's time to fit the model with a MCMC algorithm for a specific _p_ probability.

```{r}
GPDP_MCMC <- GPDPQuantReg(y ~ x, sample_data, p = 0.250)
```

Some diagnostics for the Markov Chains are available (ergodicity, autocorrelation, crosscorrelation and traces).

```{r}
diagnose(GPDP_MCMC)
```

Since it is a non-parametric model, it is focused on prediction with a credible interval. 

```{r}
predictive_data <- data_frame(x = seq(-15, 15, 0.25))
credibility <- 0.90
prediction <- predict(GPDP_MCMC, predictive_data, credibility)
```

