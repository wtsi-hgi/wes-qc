#!/usr/bin/env Rscript

list.of.packages <- c("data.table", "dplyr", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(dplyr)
library(ggplot2)

options(error = quote({
  dump.frames("errorDump", to.file=TRUE, include.GlobalEnv=TRUE); # First dump to file; this dump is not accessible by the R session.
  sink(file=stderr()); # Specify sink file to redirect all output.
  dump.frames(); # Dump again to be able to retrieve error message and write to error log; this dump is accessible by the R session since not dumped to file.
  cat(attr(last.dump,"error.message")); # Print error message to file, along with simplified stack trace.
  cat('\nTraceback:');
  cat('\n');
  traceback(2); # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
  sink()}))

read_from_json <- function(filename, type){
  data <- jsonlite::read_json(filename)
  df <- rbindlist(data[[type]])
  df[['filter']] <- names(data[[type]])
  total_tp <- data[[paste0(type,'_total_tp')]]
  total_fp <- data[[paste0(type,'_total_fp')]]
  df$TP <- 100 * df$TP / total_tp
  df$FP <- 100 * df$FP / total_fp
  df <- rename(df, precision = prec)
  return(df)
}

write("Plotting graphs to for missingness", stdout())
args <- commandArgs(TRUE)

if (length(args) != 2) {
  write("Need two arguments to run the script: <input file> <var type: snv/indel>", stdout())
  write(length(args), stdout())
  q()
}

json_data_name <- args[1]
var_type <- args[2]

#json_data_name <- 'evaluation.indel.json'
#var_type <- 'indel' #

data_basename <- sub('\\.json$', '', json_data_name)
write('Data basename', stdout())
write(data_basename, stdout())

df <- read_from_json(json_data_name, var_type)
data <- jsonlite::read_json(json_data_name)

df <- rbindlist(data[[var_type]])
df[['filter']] <- names(data[[var_type]])


df %>%
  tidyr::separate(col = 'filter', into = c('V1', 'bin', 'V2', 'DP', 'V3', 'GQ', 'V4', 'AB'), sep = '_') %>%
  select(-paste0('V', 1:4)) -> df_parsed

df %>%
  tidyr::separate(col = 'filter', into = c('V1', 'bin', 'V2', 'DP', 'V3', 'GQ', 'V4', 'AB', 'V5', 'call_rate'), sep = '_') %>%
  select(-paste0('V', 1:5)) -> df_parsed

total_tp <- data[[paste0(var_type,'_total_tp')]]
total_fp <- data[[paste0(var_type,'_total_fp')]]

df_parsed$TP <- df_parsed$TP / total_tp * 100
df_parsed$FP <- df_parsed$FP / total_fp * 100

fwrite(df_parsed, paste0(data_basename, '.csv'))

ggplot(df_parsed, aes(x=TP, y=FP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
  #xlim(c(94, 99)) + ylim(c(0, 5)) +
  facet_wrap(. ~ call_rate) + scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  ggtitle(paste0(var_type, " FP(TP)"))
ggsave(paste0(data_basename, '.FP-TP.png'))

if (var_type == 'snv') {
  ggplot(df_parsed, aes(x=t_u_ratio, y=TP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
    facet_wrap(. ~ call_rate) + scale_y_continuous(breaks = seq(0, 100, by = 5)) + scale_x_continuous(breaks = seq(0.9, 1, by = 0.005)) +
    ggtitle(paste0(var_type, " T/U ratio (TP)"))
  ggsave(paste0(data_basename,'.tu-TP.png'))
}

ggplot(df_parsed, aes(x=prec, y=recall, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
  facet_wrap(. ~ call_rate) +
  ggtitle(paste0(var_type, " recall(precision)"))
ggsave(paste0(data_basename,'.recall-prec.png'))

# indels
#ggplot(df_parsed, aes(x=TP, y=FP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
#  xlim(c(90, 95)) + ylim(c(0, 5)) + facet_wrap(. ~ call_rate)
#ggplot(df_parsed, aes(x=precision, y=recall, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
#  xlim(c(0.64, 0.81)) + ylim(c(0.64, 0.81)) + facet_wrap(. ~ call_rate)
