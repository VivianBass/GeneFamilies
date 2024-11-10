
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("exp.prof.dists_tissue calculates distances per tissue for multiple genes", {
    # Sample data with multiple genes and multiple tissues
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- c("gene1", "gene2", "gene3")
    
    # Apply the function
    distances <- exp.prof.dists_tissue(gene.accessions, expression.profiles)
    
    # Check that distances is a list with an entry per tissue
    expect_type(distances, "list")
    expect_named(distances, c("tissue1", "tissue2"))
    
    # Each distance vector should have the correct length for 3 genes (3 pairwise comparisons)
    expect_length(distances$tissue1, 3)
    expect_length(distances$tissue2, 3)
    expect_true(all(distances$tissue1 >= 0) && all(distances$tissue2 >= 0))
})

test_that("exp.prof.dists_tissue returns NA if only one gene is provided", {
    # Sample data with only one gene
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1"),
        tissue1 = c(1.2),
        tissue2 = c(2.5)
    )
    gene.accessions <- c("gene1")
    
    # Apply the function
    distances <- exp.prof.dists_tissue(gene.accessions, expression.profiles)
    
    # Expect NA as there's only one gene
    expect_true(is.na(distances))
})

test_that("exp.prof.dists_tissue handles non-matching genes gracefully", {
    # Define expression profiles with non-matching genes
    expression.profiles <- data.frame(
        FBpp_ID = c("gene4", "gene5"),
        tissue1 = c(1.0, 1.5),
        tissue2 = c(2.0, 2.5)
    )
    gene.accessions <- c("gene1", "gene2", "gene3")
    
    # Apply the function
    distances <- exp.prof.dists_tissue(gene.accessions, expression.profiles)
    
    # Expect NA as there are no matching genes
    expect_true(is.na(distances))
})

test_that("exp.prof.dists_tissue handles missing tissue columns gracefully", {
    # Sample data missing one specified tissue column
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2"),
        tissue1 = c(1.2, 2.3)
    )
    gene.accessions <- c("gene1", "gene2")
    
    # Apply the function with tissues that include a non-existent column
    expect_error(exp.prof.dists_tissue(gene.accessions, expression.profiles, tissues = c("tissue1", "tissue2")), 
                 "undefined columns selected")
})

test_that("exp.prof.dists_tissue computes distances using specified distance method", {
    # Sample data with multiple genes across tissues
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- c("gene1", "gene2", "gene3")
    
    # Apply the function with "manhattan" distance method
    distances <- exp.prof.dists_tissue(gene.accessions, expression.profiles, dist.method = "manhattan")
    
    # Check that distances vector is of correct length for each tissue
    expect_length(distances$tissue1, 3)
    expect_length(distances$tissue2, 3)
})
