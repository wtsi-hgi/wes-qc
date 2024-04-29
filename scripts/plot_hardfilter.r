#/usr/bin/Rscript

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
  sink();
  q()}))

read_from_json <- function(filename, type){
  data <- jsonlite::read_json(filename)
  df <- rbindlist(data[[type]])
  df[['filter']] <- names(data[[type]])
  df$TP <- 100 * df$TP / data$snv_total_tp
  df$FP <- 100 * df$FP / data$snv_total_fp
  df <- rename(df, precision = prec)
  return(df)
}

write("Plotting graphs to for missingness", stdout())

args <- commandArgs(TRUE)
#json_data_name <- args[1]
json_data_name <- 'evaluation.snp.json'
var_type <- 'snv' #args[2]


write("args:", stdout())
write(json_data_name, stdout())

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

df_parsed$TP <- df_parsed$TP / data$snv_total_tp * 100
df_parsed$FP <- df_parsed$FP / data$snv_total_fp * 100


ggplot(df_parsed, aes(x=TP, y=FP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
  #xlim(c(94, 99)) + ylim(c(0, 5)) + 
  facet_wrap(. ~ call_rate) + scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  ggtitle(paste0(var_type, " FP(TP)"))
ggsave(paste0(json_data_name,'_',var_type,'_FP-TP.png'))

if (var_type == 'snv') {
  ggplot(df_parsed, aes(x=t_u_ratio, y=TP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
    facet_wrap(. ~ call_rate) + scale_y_continuous(breaks = seq(0, 100, by = 5)) + scale_x_continuous(breaks = seq(0.9, 1, by = 0.005)) +
    ggtitle(paste0(var_type, " T/U ratio (TP)"))
  ggsave(paste0(json_data_name,'_',var_type,'_tu-TP.png'))
}

ggplot(df_parsed, aes(x=prec, y=recall, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
  facet_wrap(. ~ call_rate) +
  ggtitle(paste0(var_type, " recall(precision)"))
ggsave(paste0(json_data_name,'_',var_type,'_recall-prec.png'))

# indels
#ggplot(df_parsed, aes(x=TP, y=FP, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
#  xlim(c(90, 95)) + ylim(c(0, 5)) + facet_wrap(. ~ call_rate)
#ggplot(df_parsed, aes(x=precision, y=recall, color=bin, shape=DP, size=AB, alpha=GQ)) + geom_point() + theme_bw() +
#  xlim(c(0.64, 0.81)) + ylim(c(0.64, 0.81)) + facet_wrap(. ~ call_rate)
