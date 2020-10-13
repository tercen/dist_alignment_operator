library(tercen)
library(dplyr)
library(seqinr)
library(tidyr)

options("tercen.workflowId" = "a77770c3923fad0ca99b77fa8905471d")
options("tercen.stepId"     = "7db7e18c-2668-49a6-b99d-91993d0d9942")

(ctx = tercenCtx())

letters <- ctx$select(ctx$colors[[1]])[[1]]

id_na_pos <- which(is.na(ctx$cselect(ctx$cnames[[1]]))) - 1
if(length(id_na_pos) == 0) id_na_pos <- -1

df <- ctx %>% select(.ri, .ci) %>%
  mutate(letter = letters) %>%
  filter(.ci != id_na_pos) %>% 
  spread(.ci, letter)

df[is.na(df)] <- "-"

aln <- list()
aln$nam <- df[,1][[1]]
aln$seq <- as.list(apply(df[,-1], 1, function(x) paste0(x[!is.na(x)], collapse = "")))
aln$com <- NA
aln$nb <- length(aln$nam)
class(aln) <- "alignment"

sequence_type <- "protein"
if(!is.null(ctx$op.value('sequence_type'))) { sequence_type <- ctx$op.value('sequence_type') } 
matrix <- "similarity"
if(!is.null(ctx$op.value('matrix'))) { matrix <- ctx$op.value('matrix') }
gap <- 0
if(!is.null(ctx$op.value('gap'))) { gap <- as.numeric(ctx$op.value('gap')) }

if(sequence_type != "protein" & matrix == "similarity") {
  stop("Similarity matrix only available for protein sequences.")
}

dist.aln <- as.matrix(dist.alignment(aln, matrix = matrix, gap = gap))

rnames <- ctx$rselect()[[1]]

df_out <- dist.aln %>%
  as.data.frame %>%
  mutate(.ri = 1:nrow(.) - 1) %>%
  gather(.ci, dist, -.ri) %>% mutate(.ci = as.numeric(.ci))

df_out %>%
  mutate(dist_to = rnames[df_out$.ci + 1]) %>%
  ctx$addNamespace() %>%
  ctx$save()
