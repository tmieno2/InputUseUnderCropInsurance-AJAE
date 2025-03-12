# Computer programs to reproduce (though not perfectly) the results found in "Input use under crop insurance: the role of actual production history"

This repository contains the code used to generate and present the results from Mieno et al. (2018). Please note the following:

+ The code here is not the original version but a more organized and efficient rewrite. As a result, the output differs slightly from the published results. In the original version, I set a seed value in an unusual place, but I prioritized cleaner and more efficient code over exact reproducibility.

+ You still need to configure your system to use the `Rcpp` package. In hindsight, I didn't need to use C++ at all.

+ Some parallelized computations require a significant amount of RAM. Adjust the number of cores based on your system's memory capacity, or simply use `lapply()` instead.

# Reference

Mieno, T., Walters, C. G., & Fulginiti, L. E. (2018). Input use under crop insurance: The role of actual production history. American Journal of Agricultural Economics, 100(5), 1469-1485.
