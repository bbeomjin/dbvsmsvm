# dbvsmsvm
Derivative-based variable selection method for the angle-based multicategory support vector machine

```dbvsmsvm``` is an R package. ```dbvsmsvm``` provides functions to perform the derivative-based variable selection method for multicategory problems using the reinforced angle-based multicategory support vector machine (RAMSVM). It also provides functions performing the structured reinforced angle-based MSVM (SRAMSVM) and generating simulated data. 

## 1. INSTALLATION

```dbvsmsvm``` is not submitted to CRAN. Therefore, ```dbvsmsvm``` can not be installed through install.packages("dbvsmsvm") in R prompt.
Instead, ```dbvsmsvm``` can be installed through our GitHub.
Please see below to install in R.

(1) From GitHub
```{r}
> library(devtools)
> install_github("bbeomjin/dbvsmsvm")
```

## 2. USAGE NOTES

(1) Description of R functions in ```dbvsmsvm```

- Descriptions of arguments in the functions in ```dbvsmsvm``` can be obtained by help() or ? in R prompt, and documentation of ```dbvsmsvm```.   


(2) List of R functions in ```dbvsmsvm``` package

- ```dbvsmsvm``` : ```dbvsmsvm``` function is used to implement variable selection with the derivative-based variable selection method.

- ```cv.ramsvm``` : ```cv.ramsvm``` function tunes regularization parameters of the RAMSVM by k-fold cross-validation and yields the model on the complete data with optimal parameter.

- ```ramsvm``` : ```ramsvm``` is used to fit the RAMSVM with supplied hyperparameter on the given data.

- ```sramsvm``` : ```sramsvm``` function tunes hyperparameters of the SRAMSVM and then fit the SRAMSVM with optimal hyperparameters.


## 3. EXAMPLE

```{r}
# Generation of simulated data with linear decision boundary
> require(dbvsmsvm)
> n = 100; p = 10; 
> data = dbvsmsvm:::sim_gen(n = n, p = p, type = "linear")
> x = scale(data$x)
> y = data$y
> sigma = kernlab::sigest(y ~ x, scaled = FALSE)[3]


# Fit the DBVS-MSVM method with the linear and Gaussian kernel
# The number of lambda values is set to 100. 
# The number of threshold values is set to 100.
# The optimal lambda and threshold values are selected via 5-fold cross-validation with one standard error rule.
# Fit the DBVS-MSVM with the linear kernel
> dbvs_linear = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = c(2^{seq(-20, 5, length.out = 100)}),
                         Nofv = 100, kernel = "linear", scale = FALSE, cv_type = "osr", 
                         interaction = FALSE, gamma = 0.5, optModel = FALSE, nCores = 1)

# Fit the DBVS-MSVM method with the Gaussian kernel
> dbvs_radial = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = c(2^{seq(-20, 5, length.out = 100)}),
                         Nofv = 100, kernel = "gaussian", kparam = sigma, scale = FALSE, cv_type = "osr", 
                         interaction = FALSE, gamma = 0.5, optModel = FALSE, nCores = 1)

# Fit the DBVS-MSVM with the Gaussian kernel for selecting second-order interaction
dbvs_interaction = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = c(2^{seq(-20, 5, length.out = 100)}),
                       Nofv = 100, kernel = "gaussian", kparam = sigma, criterion = "0-1", scale = FALSE,
                       cv_type = "osr", interaction = TRUE, gamma = 0.5, optModel = TRUE, nCores = 1)


# Fit the SRAMSVM with the linear kernel
> sram_linear = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                        lambda_seq = 2^{seq(-20, 5, length.out = 100)},
                        lambda_theta_seq = 2^{seq(-20, 5, length.out = 100)},
                        kernel = "linear", scale = FALSE, criterion = "0-1",
                        isCombined = TRUE, cv_type = "osr", nCores = 1)

# Fit the SRAMSVM with the Gaussian kernel
> sram_radial = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                        lambda_seq = 2^{seq(-20, 5, length.out = 100)},
                        lambda_theta_seq = 2^{seq(-20, 5, length.out = 100)},
                        kernel = "gaussian", kparam = sigma, scale = FALSE, criterion = "0-1",
                        isCombined = TRUE, cv_type = "osr", nCores = 1)
                      
# Fit the SRAMSVM with the Gaussian kernel with second-order interaction
> sram_radial_interaction = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                                   lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                                   lambda_theta_seq = 2^{seq(-20, 0, length.out = 100)},
                                   kernel = "gaussian", kparam = sigma, scale = FALSE, criterion = "0-1",
                                   isCombined = TRUE, cv_type = "osr", nCores = 1)

```

