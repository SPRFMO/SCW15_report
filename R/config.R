# config.R - DESC
# jmMSE/config.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# LOAD packages

library(TAF)
library(mse)
library(FLjjm)
library(qs)

# SOURCE extra functions
source(here::here("R", "utilities.R"))
source(here::here("R", "plots.R"))

# SETUP parallel
library(doFuture)

options(future.globals.maxSize=1500 * 1024 ^ 2)

# SET progressr reporting
#handlers(global=TRUE)

# --- OPTIONS

# SELECT tuning period
tperiod <- seq(2034, 2042)

# SET MSE arguments, initial year, frequency, final year
mseargs <- list(iy=2024, frq=1, fy=2045)

# LOAD performance statistics
data(statistics, package="FLjjm")

# SET metrics
mets <- list(C=function(x) catch(x), F=function(x) fbar(x),
  SB=function(x) ssb(x))
