# TIMSS 2023 grade 4 analysis of the paper with the package interface.
# Input: the file written by timss/prepare_timss2023.R of the replication
# materials (responses, SENWGT weights, countries, booklets, item metadata).
#   Rscript timss2023.R output/timss/timss2023_G4_2domain.rds
library(subsampleIRT)
args <- commandArgs(trailingOnly = TRUE)
tim <- readRDS(args[1])

dat <- irt_data(tim$Y, group = tim$country[tim$group],
                booklet = sub("^\\S+ ", "", tim$block), weights = tim$weight)
mod <- irt_model(dat, dimension = tim$dimension, guessing = tim$free_guessing,
                 means = "sum_zero", nodes = 11)
print(dat); print(mod)

ctrl <- list(chunk = 5000L)
des <- sub_design(dat, mod, targets = c("difficulty", "guessing", "correlation", "mean", "joint"),
                  pilot_size = 10000, seed = 20261400, control = ctrl)
print(des)

# one target-specific estimate of the item difficulties with standard errors
fd <- sub_fit(des, size = 30000, target = "difficulty", hessian = "forward", seed = 1)
print(fd)
print(head(summary(fd)))

# the joint design, with intervals for all item parameters and country means
fj <- sub_fit(des, size = 30000, target = "joint", hessian = "forward", seed = 2)
summary(fj, target = "mean")

# guessing-penalty path on one joint-target draw
path <- sub_path(des, size = 30000, target = "joint", seed = 3,
                 penalties = lapply(c(0.002, 0.005, 0.010, 0.020), function(s) irt_penalty(guessing = s)))
sapply(path, function(f) mean(target_estimate(f, "guessing")))

# full-data fit for comparison (about 10 minutes on four cores)
# full <- irt_fit(dat, mod, start = des$pilot$eta, control = ctrl)
