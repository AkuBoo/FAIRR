# Description
Functional Approximation of Impulse Responses (FAIR) on R. 

Follows [Barnichon and Matthes (2018) *Functional Approximation of Impulse Responses*,](https://www.sciencedirect.com/science/article/abs/pii/S0304393218302356) JME., using Gaussian basis functions. The functions are estimated using [Stan](https://mc-stan.org/)'s No-U-Turn-Sampler (NUTS) algorithm.

This code is **NOT** an exact replica of the original MATLAB code by C. Matthes, available [here](https://github.com/cm1518/FAIR). I am not the best coder out there so use at your own risk, check the code if unconvinced. [Link to the article on R. Barnichon's drive](https://drive.google.com/file/d/10wjsWyvBQHSryPQAA0VO8Tc5DCVnaE63/view?pli=1). PLEASE notify me if you see issues. 

Setting up the function is simple, install the 3 files in your working directory and put `source('FAIRR.R')` in the console. The function `fair()` can be used after. 

*Required* :
- A `R` installation (obviously). 
- `install.packages(c("vars", "rstan", "pracma", "parallel", "R.utils"))`
- A C++ installation for `rstan`.  
## Usage
`fair(Y, Z, X = NULL, H, K = 1, asym = FALSE, maxvarlag = 12, varci = FALSE, varboot = 500, varout = FALSE, spanloess = 0.5, sdpar = 10, cores = NA, chains = 4, iter = 2000, adapt_delta = 0.95, max_treedepth = 12, timeout = 180)`

Only Y, Z and H are mandatory arguments.
## Arguments
`Y, Z, X` get converted to matrices, they need matching lengths (number of rows).
- `Y` Dependent variable (the vector we want the IRF of).
- `Z` Shocks (the vector we want the IRF wrt to).
- `X` controls.
- `H` Horizon for the IRF
- `K` Number of Gaussian basis functions, currently only `1` or `2` can be specified there. `1` by default : 
$$\psi(h) = \sum^K_{k=1} a_k \exp\left(- \frac{(h-b_k)^2}{c_k^2} \right) \qquad \forall h \geq 0$$
- `asym` Should there be distinct IRF's for positive and negative orthogonal shocks ? `FALSE` by default.
- `maxvarlag` Maximum number of lags for the VAR($p$) model. Can be reduced to improve performance.
- `varci` Should the VAR IRF's confidence intervals be computed ? *Note : `Error in obs - p : non-numeric argument to binary operator`* *is a common error if there are too few observations.*
- `varboot` Bootstrap parameter for `vars::irf()`. Only applicable when `varci = TRUE`.
- `varout` Should the function output the `vars::VAR()` object ?
- `spanloess` See `span` in the `loess()` function. Controls how smooth the VAR's IRF is made before extracting priors.
- `sdpar` Standard deviation of the parameters `b` and `c` when estimating using `stan` (the parameters are assumed to be normally distributed around the prior with standard deviation `sdpar`). in the `.stan` files. 
- `cores` Amount of CPU cores used when estimating. Use at your own risk on laptops.
- `chains, iter, adapt_delta, max_treedepth` : see `stan()` documentation.
- `timeout` Maximum time spent on **each** `stan` estimation (if `asym = TRUE` it doubles!). 
## Value
A `list` :
	`$fair` FAIR output
		**If `asym = FALSE`**
		`$params` parameters
			`$mean` mean
			`$low` 2.5 centile
			`$high` 97.5 centile
		`$psi` posterior distribution
			`$mean` mean
			`$low` 2.5 centile
			`$high` 97.5 centile
		**If `asym = TRUE`**
		`$positive`
			`$params` parameters
				`$mean` mean
				`$low` 2.5 centile
				`$high` 97.5 centile
			`$psi` posterior distribution
				`$mean` mean
				`$low` 2.5 centile
				`$high` 97.5 centile
		`$negative`
			`$params` parameters
				`$mean` mean
				`$low` 2.5 centile
				`$high` 97.5 centile
			`$psi` posterior distribution
				`$mean` mean
				`$low` 2.5 centile
				`$high` 97.5 centile
	`$var` VAR estimated 
		`$irf` Impulse-Response Function
			`$mean` mean
			`$low` 2.5 centile
			`$high` 97.5 centile
		`$shocks` Orthogonalized shocks
			`$all` Matrix of all the shocks
		`$priors` Priors extracted from the smoothed IRF
			**If K=1**
			`$a, $b, $c` with `$mean, $low, $high` 
			**If K=2**
			`$a1, $b1, $c1, $a2, $b2, $c2` with `$mean, $low, $high`
	`$stan` The `stan::summary` output with default options, see their documentation.
	`$data` Data used
		`$Y`
		`$Z`
		`$X`
## Known issues
If the noise to signal ratio is high, the output can and probably will be garbage. I've tried to compare VAR IRF's with FAIR on simulated processes and FAIR seemed consistently better than VAR in highly noisy settings but Stan threw a lot of warnings as the priors for the FAIR come from the VAR's IRF and in these settings the posterior distribution can be quite different from the prior distribution. Estimation can be extremely slow if the noise to signal ratio is high. 
## Files
- `FAIRR.R`. Contains the code for the function `fair()`.
- `FAIRone.stan` Contains the `stan` code for the estimation of one Gaussian.
- `FAIRtwo.stan` Contains the `stan` code for the estimation of two Gaussians.
## Planned future updates 
More basis functions, support for 3 and 4 Gaussians. Keep in mind interpretability matters so using more than 2 Gaussians starts to defeat the point. I will probably add more options for the `stan()` function. 
