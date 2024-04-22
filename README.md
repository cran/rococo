# RoCoCo - An R Package Implementing a Robust Rank Correlation Coefficient and a Corresponding Test

Classical rank correlation measures are designed for ordinal data, but they are not ideally suited for measuring rank correlation for numerical data that are perturbed by noise. We introduced a family of robust rank correlation measures on the basis of Goodman's and Kruskal's gamma, where the classical ordering of real numbers is replaced by some fuzzy ordering with smooth transitions - thereby ensuring that the correlation measure is continuous with respect to the data. The present R package implements this rank correlation measure along with a corresponding statistical test. Moreover, the package also implements the Gaussian rank correlation estimator and a corresponding statistical test.

Although the package is maintained by Ulrich Bodenhofer, large parts of the
package itself have been implemented by Martin Krone.

## Installation

The package can be installed from
[CRAN](https://CRAN.R-project.org/package=rococo). Therefore, the the simplest way to install the package is to enter
```
install.packages("rococo")
```
into your R session. If, for what reason ever, you prefer to install the package manually, follow the instructions in the [user manual](https://cran.r-project.org/package=rococo/vignettes/rococo.pdf).

## User support

If you encounter any issues or if you have any question that might be of interest also for other users, before writing a private message to the package developers/maintainers, please create an issue in this repository and also consider posting to the [R-help Mailing List](https://stat.ethz.ch/mailman/listinfo/r-help) or on [StackOverflow](https://stackoverflow.com/). For other matters regarding the package, please contact the package author.

## Citing this package

If you use this package for research that is published later, you are kindly asked to cite it as follows:

- U. Bodenhofer, M. Krone, and F. Klawonn (2013). Testing noisy numerical data for monotonic association. *Inform. Sci.* **245**:21-37. DOI: [10.1016/j.ins.2012.11.026](http://doi.org/10.1016/j.ins.2012.11.026).
- U. Bodenhofer and F. Klawonn (2008). Robust rank correlation coefficients on the basis of fuzzy orderings: initial steps. *Mathware Soft Comput.* **15**(1):5-20.
