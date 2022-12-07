# library(oiff); library(testthat); source("test-findOptimalFilter.R")

reference <- function(pvalues, covariates, threshold=0.05) {
    best.hits <- -1L
    best.filter <- 0

    for (x in covariates) {
        keep <- covariates <= x
        curp <- pvalues[keep]
        nhits <- sum(p.adjust(curp, method="BH") <= threshold)
        if (nhits > best.hits || (nhits == best.hits && x > best.filter)) {
            best.hits <- nhits
            best.filter <- x
        }
    }

    list(threshold = best.filter, number = best.hits)
}

test_that("findOptimalFilter works when compared to a reference", {
    set.seed(999)
    pvalues <- c(runif(900), rbeta(100, 1, 50))
    filter <- c(rnorm(900), rnorm(100) + 2)

    ref <- reference(pvalues, -filter)
    ref$threshold <- -ref$threshold
    expected <- findOptimalFilter(pvalues, filter, above=TRUE)

    expect_identical(ref, expected)
})

test_that("findOptimalFilter is pretty fast", {
    set.seed(1111)
    pvalues <- c(runif(999000), rbeta(1000, 1, 50))
    filter <- c(rnorm(999000), rnorm(1000) + 2)
    expected <- findOptimalFilter(pvalues, filter, above=TRUE)
    expect_type(expected, "list")
})
