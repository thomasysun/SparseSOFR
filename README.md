# SparseSOFR

`SparseSOFR` implements a Bayesian scalar-on-function regression model with a sparsely observed and/or noisy functional covariate.

## Background

Functional and longitudinal data are commonly measured with error, with the latter frequently encountering issues of missingness and sparse measurements. We consider the following scalar-on-function regression model along with a measurement error model for the functional covariate.

$$
y_i = \alpha_0 + \boldsymbol{w}_i^{\prime}\boldsymbol{\alpha} + \int_{T} X_i(\tau)\phi(\tau)d\tau +\epsilon_{yi} ,  \quad  \epsilon_{yi}(\tau) \overset{iid}{\sim} N(0, \sigma^2_{\epsilon_y})
$$

$$
x_i(\tau) = X_i(\tau)  + \epsilon_{xi}(\tau) \\ 
= \sum^{K}_{k=1}f_k(\tau_{})\beta_{k,i} + \epsilon_{xi}(\tau) ,  \quad  \epsilon_{xi}(\tau) \overset{iid}{\sim} N(0, \sigma^2_{\epsilon_x})
$$

$y_i$ is a scalar response variable, $\boldsymbol{w}_i$ contains scalar covariates with intercept $\alpha_0$ and coefficients $\boldsymbol{\alpha} = \{\alpha_1, \dots, \alpha_L\}^{\prime}$, and $X_i(\tau)$ is the functional covariate with functional coefficient $\phi(\tau)$. The coefficient $\phi(\tau)$ explains the effect of the functional predictor at $\tau$ on the response value. 

In practice, we do not get to observe the functional covariate $X_i(\tau)$, rather we observe the discrete vector of noisy measurements of a function, $x_i$. We use a functional mixed effects model to pool sparse observations for curve fitting. A functional factor model is used to model the unknown basis functions, which allows for estimating high sparse functional observations in a data-driven fashion and incorporates uncertainty about the unknown basis.


## Installation

You can install and load the package using the following lines:
```{r install}
# install.packages("devtools")
# devtools::install_github("thomasysun/SparseSOFR")
library(SparseSOFR) 
```

## Usage

The primary function `ssofr()` runs the scalar-on-function regression model. It requires three inputs at minimum:

- `y` : vector of scalar responses.
- `x` : matrix of the functional responses. Each column corresponds to one function and each row is a measurement at a single time point on a common grid. May contain NAs.
- `Tau` : vector specifying common grid of time points on which function is observed.

Scalar covariates may also be added via a design matrix `w`. See documentation for full details.

Sparse and/or noisy curves may also be fit in a standalone model using `sffm()`.



