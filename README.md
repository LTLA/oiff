# Optimizing an independent filter for FDR control

## Overview

Given a set of p-values and a filter statistic that is independent of the p-values under the null,
the **oiff** library identifies the filter threshold that maximizes the number of discoveries at a given FDR threshold.
Conceptually, it yields the same result as the following naive procedure:

1. Retain only those hypotheses where the filter statistic is below some filter threshold.
2. Apply the Benjamini-Hochberg (BH) method to the retained hypotheses.
3. Count the number of discoveries among the retained hypotheses at a given FDR threshold.
4. Repeat 1-3 to find the filter threshold that maximizes the number of discoveries.

This can provide a "sensible" choice for the filter threshold when no _a priori_ setting is available.
For example, we often filter out low-abundance features prior to differential analyses of genomic data,
on the basis that the abundance of a feature is usually independent of its p-value.

## Quick start

C++ users can just link to [the header](include/oiff/oiff.hpp) and run:

```cpp
#include "oiff/oiff.hpp"

std::vector<double> pvalues; // fill with p-values
std::vector<double> covariates; // fill with covariates

// Finds the optimal filter at a FDR threshold of 0.05.
oiff::OptimizeFilter runner;
runner.fdr_threshold = 0.05;
auto res = runner.run(pvalues.size(), pvalues.data(), covariates.data());
res.middle; // one choice of filter threshold
res.number; // number of discoveries

// Run with subsampling and take the average of subsample iterations.
auto res2 = runner.run_subsample(pvalues.size(), pvalues.data(), covariates.data());
double mean_threshold = 0;
for (auto x : res2) {
    mean_threshold += x.middle;
}
mean_threshold /= res2.size();
```

R users can install [the test package](R/) and run the example:

```r
library(oiff)
pvalues <- c(runif(9900), rbeta(100, 1, 50))
filter <- c(rnorm(9900), rnorm(100) - 2)
findOptimalFilter(pvalues, filter)
```

Check out the [reference documentation](https://ltla.github.io/oiff) for more details.

## Building projects 

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```
include(FetchContent)

FetchContent_Declare(
  oiff
  GIT_REPOSITORY https://github.com/LTLA/oiff
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(oiff)
```

Then you can link to **oiff** to make the headers available during compilation:

```
# For executables:
target_link_libraries(myexe oiff)

# For libaries
target_link_libraries(mylib INTERFACE oiff)
```

## Comments on performance

Computationally, **oiff** uses an interval tree to avoid repeated invocations of the BH method.
This means that the algorithm is very fast for large numbers of hypotheses:

```r
library(oiff)
pvalues <- c(runif(999000), rbeta(1000, 1, 50))
filter <- c(rnorm(999000), rnorm(1000) + 2)
system.time(expected <- findOptimalFilter(pvalues, filter, above=TRUE))
##    user  system elapsed
##   0.458   0.008   0.466
```

Statistically, this approach is flawed as it does not guarantee control of the FDR.
By allowing the filter threshold to vary in a manner that depends on the p-values,
**oiff** will systematically include more false discoveries than allowed for under the BH method.
Here is a simple demonstration of the problem:

```r
library(oiff)
num.discoveries <- numeric(1000)
ref.discoveries <- numeric(1000)

for (it in seq_along(num.discoveries)) {
    # Generating null hypotheses.
    pval <- runif(100) 
    filter <- rnorm(100)

    # Injecting a single true positive that is always retained.
    pval <- c(0, pval)
    filter <- c(100, filter)

    # Using an optimal filter threshold.
    expected <- findOptimalFilter(pval, filter, threshold=0.05, above=TRUE)
    num.discoveries[it] <- expected$number

    # Compared to a constant filter.
    above.zero <- pval[filter >= 0]
    ref.discoveries[it] <- sum(p.adjust(above.zero, method="BH") <= 0.05)
}

# Calculating the FDR after removing the lone true positive: 
mean((num.discoveries - 1) / num.discoveries)
## [1] 0.1866667
mean((ref.discoveries - 1) / ref.discoveries)
## [1] 0.04925
```

A practical mitigation is to derive the threshold from a small subsample of hypotheses.
This preserves any dependencies between the p-values and filter statistic _under the alternative hypothesis_,
thus ensuring that we still reap the benefits of filter optimization.
The use of a small subsample means that the chosen filter threshold is independent of the p-values for the remaining hypotheses,
limiting the severity of the loss of FDR control (assuming that the various hypotheses are independent of each other).
This is inspired by the cross-validation procedure in the [**IHW**](https://bioconductor.org/packages/IHW) package.

```r
library(oiff)
num.discoveries <- numeric(1000)

for (it in seq_along(num.discoveries)) {
    # Generating null hypotheses.
    pval <- runif(100) 
    filter <- rnorm(100)

    # Injecting a single true positive that is always retained.
    pval <- c(0, pval)
    filter <- c(100, filter)

    # Using an optimal filter threshold based on a subsample. 
    expected <- findOptimalFilter(pval, filter, threshold=0.05, above=TRUE, subsample=0.1)
    keep <- filter >= expected$middle
    num.discoveries[it] <- sum(p.adjust(pval[keep], method="BH") <= 0.05)
}

mean((num.discoveries - 1) / num.discoveries)
## [1] 0.05008333
```
