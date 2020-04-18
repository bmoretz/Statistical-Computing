library(data.table)
library(RCurl)
library(XML)
library(stringr)
library(tidyr)
library(here)

data.dir <- file.path(here::here(), "Case Studies", "datasets")

Url <- "http://www.cherryblossom.org/results/1999/cb99f.html"


theurl <- Url

webpage <- getURL(theurl)
webpage <- readLines(tc <- textConnection(webpage)); close(tc)

pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)

pre_data <- pagetree["//pre"]

lines <- sapply(pre_data, xmlValue)

table_long <- strsplit(lines, "\n")[[1]]

total_rows <- length(table_long)
row_data <- data.table(values = table_long[-1], stringsAsFactors = F)

num_rows <- nrow(row_data)
to_parse <- row_data[3:num_rows]

parse_row <- function(row) {
  
  row_num <- as.numeric(substr(row, 1, 5))
  
  place <- substr(row, 6, 15)
  full_name <- substr(row, 16, 38)
  age <- substr(row, 38, 40)
  hometown <- substr(row, 41, 60)
  time <- substr(row, 61, 68)
  pace <- substr(row, 69, 76)
  
  return(data.frame(row_num, place, full_name, age, hometown, time, pace))
}

parsed <- apply(to_parse, 2, parse_row)

results <- data.table::rbindlist(parsed)


test_row <- row_data[3, 1]
test_row

substr(test_row, 1, 5)
substr(test_row, 6, 15)
substr(test_row, 16, 38)
substr(test_row, 38, 40)
substr(test_row, 41, 60)
substr(test_row, 61, 68)
substr(test_row, 69, 76)
