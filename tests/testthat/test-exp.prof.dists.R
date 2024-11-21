
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("exp.prof.dists calculates Euclidean distances for multiple genes", {
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        Species = c("spec1", "spec1", "spec1"),  
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    
    expect_length(unlist(distances), 3)  
    expect_true(all(unlist(distances) >= 0))
})

test_that("exp.prof.dists returns NA if only one gene is provided", {
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1"),
        Species = c("spec1"),  
        tissue1 = c(1.2),
        tissue2 = c(2.5)
    )
    gene.accessions <- list("gene1")
    
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    expect_true(all(is.na(unlist(distances))))  
})

test_that("exp.prof.dists handles missing columns gracefully", {
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2"),
        Species = c("spec1", "spec1"),  
        tissue1 = c(1.2, 2.3)
    )
    gene.accessions <- list("gene1", "gene2")
    
    expect_error(exp.prof.dists(gene.accessions, expression.profiles, 
                               tissues = c("tissue1", "tissue2")))
})

test_that("exp.prof.dists computes distances using specified distance method", {
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        Species = c("spec1", "spec1", "spec1"),  
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    distances <- exp.prof.dists(gene.accessions, expression.profiles, dist.method = "manhattan")
    expect_length(unlist(distances), 3)  
})

test_that("exp.prof.dists returns empty list if no genes match", {
    expression.profiles <- data.frame(
        FBpp_ID = c("gene4", "gene5"),
        Species = c("spec1", "spec1"),  
        tissue1 = c(1.0, 1.5),
        tissue2 = c(2.0, 2.5)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    expect_equal(length(distances), 0)  
})
