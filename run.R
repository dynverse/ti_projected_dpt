#!/usr/local/bin/Rscript

task <- dyncli::main()

# load libraries
library(dyncli, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

library(destiny, warn.conflicts = FALSE)

#####################################
###           LOAD DATA           ###
#####################################
expression <- task$expression %>% as.matrix
params <- task$params
priors <- task$priors

start_cell <-
  if (!is.null(priors$start_id)) {
    sample(priors$start_id, 1)
  } else {
    NULL
  }

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################
# run diffusion maps
dm <- destiny::DiffusionMap(
  data = expression,
  sigma = params$sigma,
  distance = params$distance,
  n_eigs = params$ndim,
  density_norm = params$density_norm,
  n_local = params$n_local,
  vars = params$features_id
)

# run DPT
if (!is.null(start_cell)) {
  tips <- which(rownames(expression) %in% start_cell)
} else {
  tips <- destiny::random_root(dm)
}
dpt <- destiny::DPT(
  dm,
  w_width = params$w_width,
  tips = tips
)

# find DPT tips
tips <- destiny::tips(dpt)
tip_names <- rownames(expression)[tips]

# TIMING: done with trajectory inference
timings$method_aftermethod <- Sys.time()

#   ____________________________________________________________________________
#   Save output                                                             ####

cell_ids <- rownames(expression)

# construct grouping
grouping <- dpt@branch[,1] %>%
  ifelse(is.na(.), 0, .) %>%
  as.character() %>%
  paste0("Tip", .) %>%
  set_names(cell_ids)

# retrieve dimred
dimred <- dpt@dm@eigenvectors %>%
  magrittr::set_colnames(., paste0("Comp", seq_len(ncol(.)))) %>%
  magrittr::set_rownames(cell_ids)

# get cluster assignment
grouping <- dpt@branch[,1] %>%
  ifelse(is.na(.), 0, .) %>%
  as.character()
branches <- sort(unique(grouping))

# calculate cluster medians
dimred_milestones <- t(sapply(branches, function(br) colMeans(dimred[grouping == br,,drop=F])))

# create star network
milestone_network <- tibble(
  from = "0",
  to = setdiff(branches, "0"),
  length = sqrt(rowMeans((dimred_milestones[from,] - dimred_milestones[to,])^2)),
  directed = FALSE
)

# save output
output <-
  wrap_data(cell_ids = cell_ids) %>%
  add_dimred_projection(
    dimred = dimred,
    grouping = grouping,
    dimred_milestones = dimred_milestones,
    milestone_network = milestone_network
  ) %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
