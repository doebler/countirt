# countirt
 R package for count data IRT models.

This package implements methods for count data Item Response Theory (IRT), specifically:

- the Two-Parameter Conway-Maxwell-Poisson model (2PCMPM; Beisemann, 2022)
- the Two-Parameter Poisson Counts Model (2PPCM; Myszkowski & Storme, 2021)
- the Distributional Regression Test Model (DRTM) and the Count Latent Regression Model (CLRM) which are explanatory count IRT models based on the 2PCMPM (Beisemann, Forthmann, & Doebler, 2024)
- the Poisson variants of the DRTM and the CLRM

This package is currently in early development and will be continuously extended. At this development stage, there might be bugs. Checks and warnings so far are only minimal, requiring users to be mindful themselves of whether their data and results seem plausible. Please be aware of that if you use the package. So far, you can fit all implemented models (including some constrained versions of these models) with the cirt function. Please consult the documentation of the cirt function. 

Please be aware that the CMP models and algorithms were developed for and assessed in the data settings in Beisemann (2022) and Beisemann et al. (2022). If you wish to use them in different data settings (e.g., more strongly skewed responses, larger counts, more extreme parameter values, such as larger discriminations / slopes, etc.), assessing their performance first in simulations is strongly recommended. In the 2PCMPM / DRTM / CLRM algorithms as implemented in countirt, CMP quantities are obtained via interpolation from a look up grid (set up for mus between 0.001 and 200, and for nus between 0.01 and 50). Parameter estimation accuracy will likely be best if in all iterations of the algorithm for all parameters, the CMP quantities that need to be interpolated stay well within the bounds of this table. The more extreme the estimates are, the further one will move towards the borders of the look up table, in all likelihood decreasing parameter estimation accuracy. Similarly, as the algorithm was developed for settings with overall smaller counts (see Beisemann, 2022, and Beisemann et al., 2022), the grid truncates expected counts (mus) at 200 before interpolation (with no extrapolation carried out). Larger counts than that are thus not going to be accounted for accurately by the algorithm. 

You can install the package from GitHub e.g. using the devtools package with devtools::install_github("mbsmn/countirt"). Please note that the package includes C++ code which is tied into R using the Rcpp and RcppGSL packages. In C++, I use the GSL library. In order to be able to install and use the countirt package smoothly, you need to install GSL on your machine (this is not an R package, but a C/C++ library). If you have a Mac, you can do so e.g. with homebrew. There are tutorials online of how you can install the GSL library. You only need to have it, the rest should be taken care of by Rcpp and RcppGSL.

In the future, I hope to provide vignettes and some examples here for how to use the countirt package. Until then, you can have a look at the example script here (https://osf.io/dzcyt/) for some code using the countirt package.

If you have any questions, you can reach me at: beisemann@statistik.tu-dortmund.de

References:

Beisemann, M. (2022). A flexible approach to modeling over-, under-and equidispersed count data in IRT: The two-parameter Conway-Maxwell-Poisson model. British Journal of Mathematical and Statistical Psychology, (Advanced online publication). https://doi.org/10.1111/bmsp.12273

Beisemann, M., Forthmann, B., & Doebler, P. (2022). Understanding Ability and Reliability Differences Measured with Count Items: The Distributional Regression Test Model and the Count Latent Regression Model. Multivariate Behavioral Research, (Advance Online Publication), 1–21. https://doi.org/10.1080/00273171.2023.22885

Myszkowski, N., & Storme, M. (2021). Accounting for variable task discrimination in divergent thinking fluency measurement: An example of the benefits of a 2-Parameter Poisson Counts Model and its bifactor extension over the Rasch Poisson Counts Model. The Journal of Creative Behavior, 55 (3), 800–818. https://onlinelibrary.wiley.com/doi/abs/10.1002/jocb.490

