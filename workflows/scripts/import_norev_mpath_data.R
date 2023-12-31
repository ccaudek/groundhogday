# Script name: import_norev_mpath_data.R
# Project: groundhogday
# Script purpose: For sessions without a reversal, combine the all the 
# individual subjects' Excel files into a single file.
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed Jun  7 22:16:48 2023
# Last Modified Date: Fri Oct 13 09:20:41 2023
#
# 👉 Note that it is necessary to refer to "Task2noRev_1".

suppressPackageStartupMessages(library("rio"))
suppressPackageStartupMessages(library("tidyverse"))

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


# Read subject codes from Excel file -------------------------------------------

print("Loading data")
subj_codes <- rio::import(
  # snakemake@input[["xlsx"]]
  "data/raw/data.xlsx"
) |>
  dplyr::select(
    subj_code_1, Task2noRev_1,
    context_1, control_1,
    post_context_1, coin, TIME_start,
    TIME_end, TIME_total
  )

subj_codes <- subj_codes[!is.na(subj_codes$Task2noRev_1), ]
subj_codes$subj_idx <- subj_codes$Task2noRev_1
subj_codes$Task2noRev_1 <- NULL


# Import no-reversal PRL data form multiple files ------------------------------

# dir <- snakemake@input[["dir_data"]]
dir <- "data/raw/experiment_data"

# Filenames for sessions without reversal:
file_names <-
  as.character(list.files(path = dir, pattern = "prl_EMA_noreversal"))

n_files <- length(file_names)
n_files
d_list <- list()

for (i in 1:n_files) {
  d <- read.table(stringr::str_glue(dir, "/", file_names[i]))

  d$subj_idx <- file_names[i]
  d$epoch <- d$V1
  d$item_number <- d$V3
  d$rt <- d$V4
  d$has_answered <- d$V5
  d$image_more_rewarded <- d$V6
  d$rt_choice <- d$V7
  d$feedback <- ifelse(d$V8 == 1, 1, 0)
  d$image_chosen <- d$V9
  # 1: mi comporterei nello stesso modo
  # 2: mi comporterei in modo diverso
  d$intertrial_delay <- d$V10
  d$img_1_position <- d$V11 # 230 = dx; -230 sn
  d$img_2_position <- d$V12 # 230 = dx; -230 sn
  d$total_feedback <- d$V13 # mood from the starting value + feedbacks

  d$trial <- 1:30

  d_list[[i]] <- d
}

# Convert list into data.frame
df <- do.call(rbind.data.frame, d_list)


# Join PRL data and personal information ---------------------------------------

d <- left_join(df, subj_codes, by = "subj_idx")


# Save for later use -----------------------------------------------------------

print("Saving output")
# saveRDS(d, file = snakemake@output[["rds"]])
saveRDS(d, file = "data/prep/groundhog_norev_raw.RDS")

# eof ----
