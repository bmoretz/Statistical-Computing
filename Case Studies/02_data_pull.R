library(data.table)
library(RCurl)
library(XML)
library(stringr)
library(tidyr)
library(here)
library(urltools)

data.dir <- file.path(here::here(), "Case Studies", "datasets")

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

url <- "http://www.cherryblossom.org/results/2003/CB03-M.HTM"
rm(url)

save_data <- function(url) {
  
  url <- tolower(url)
  parts <- unlist(strsplit(url, "/"))
  year <- parts[5]
  
  mf <- str_detect(unlist(strsplit(parts[6], "\\."))[1], "m")
  
  file_path <- data.dir
  
  if(mf) {
    file_path <- file.path(file_path, "race_men", paste0(year, "m", ".csv"))
  } else {
    file_path <- file.path(file_path, "race_women", paste0(year, "f", ".csv"))
  }

  webpage <- getURL(url)
  webpage <- readLines(tc <- textConnection(webpage)); close(tc)
  
  pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)
  
  pre_data <- pagetree["//pre"]
  
  lines <- sapply(pre_data, xmlValue)
  
  table_long <- strsplit(lines, "\n")[[1]]
  
  total_rows <- length(table_long)
  row_data <- data.table(values = table_long[-1], stringsAsFactors = F)
  
  num_rows <- nrow(row_data)
  to_parse <- row_data[3:num_rows]
  
  parsed <- apply(to_parse, 2, parse_row)
  
  results <- data.table::rbindlist(parsed)
  
  data.table::fwrite(file = file_path, results)
}

ubase <- "http://www.cherryblossom.org/"
url <- paste(ubase, "results/1999/cb99m.html", sep = "")

extractResTable <- function(url, year ) {
  
  print(paste("processing:", url, "race results for year:", year, ","))
  
  doc <- htmlParse(url, encoding="UTF-8")
  
  if (year == 2000) {
    # Get text from 4th font element
    # File is ill-formed so <pre> search doesn't work
    ff <- getNodeSet(doc, "//font")
    txt <- xmlValue(ff[[4]])
  } else {
    txt <- xpathSApply(doc, '//pre', xmlValue)
  }
 
  result <- character(length = 0L)
  
  if(grepl("\n", txt, fixed = TRUE) == T) {
    result <- strsplit(txt, "\\r\\n")[[1]]
  } else if(grepl("\n", txt, fixed = TRUE) == T) {
    result <- strsplit(txt, "\\n")[[1]]
  } else {
    result <- txt
  }
  
  return(result)
}

mens_result_urls <-c(
  "results/1999/cb99m.html",
  "results/2000/Cb003m.htm",
  "results/2001/oof_m.html",
  "results/2002/oofm.htm",
  "results/2003/CB03-M.HTM",
  "results/2004/men.htm",
  "results/2005/CB05-M.htm",
  "results/2006/men.htm",
  "results/2007/men.htm",
  "results/2008/men.htm",
  "results/2009/09cucb-M.htm",
  "results/2010/2010cucb10m-m.htm",
  "results/2011/2011cucb10m-m.htm",
  "results/2012/2012cucb10m-m.htm")

womens_result_urls <- c(
  "results/1999/cb99m.html",
  "results/2000/Cb003f.htm",
  "results/2001/oof_f.html",
  "results/2002/ooff.htm",
  "results/2003/CB03-F.HTM",
  "results/2004/women.htm",
  "results/2005/CB05-F.htm",
  "results/2006/women.htm",
  "results/2007/women.htm",
  "results/2008/women.htm",
  "results/2009/09cucb-F.htm",
  "results/2010/2010cucb10m-f.htm",
  "results/2011/2011cucb10m-f.htm",
  "results/2012/2012cucb10m-f.htm"
)

save_txt_data <- function(url, year, class) {
  data <- extractResTable(url, year )
  
  file_path <- file.path(data.dir, 
                         paste(class, "txt", sep = "_"),
                         paste(year, ".txt", sep = ""))
  
  result <- tryCatch({
    
    file_conn <- file(file_path)
    writeLines(data, file_conn)
    close(file_conn)
    
    return(TRUE)
  }, error = function(e){
    return(e)
  })
}

save_class_result <- function( info ) {
  
  class <- head(info$class, 1)
  
  result_tables <- mapply(extractResTable, 
                          url = info$url, 
                          year = info$year)
  
  names(result_tables) <- info$year
  
  result_path <- file.path(data.dir, paste(class, ".RDS", sep = ""))
  
  saveRDS(result_tables, file = result_path)
}

race_info_women <- data.frame(url = paste(ubase, mens_result_urls, sep=""), year = 1999:2012, class = rep("women", 14))
race_info_men <- data.frame(url = paste(ubase, mens_result_urls, sep=""), year = 1999:2012, class = rep("men", 14))

race_info <- rbind(race_info_men, race_info_women)

save_class_result(race_info_men)

sapply(result_tables, length)

# save all race results out as plain text files.

mapply(save_txt_data, 
       url = race_info$url, 
       year = race_info$year, 
       class = race_info$class)


# save_data("http://www.cherryblossom.org/results/1999/cb99f.html")
# save_data("http://www.cherryblossom.org/results/1999/cb99m.html")

# save_data("http://www.cherryblossom.org/results/2000/Cb003m.html")
# save_data("http://www.cherryblossom.org/results/2000/Cb003f.html")

#save_data("http://www.cherryblossom.org/results/2001/oof_m.html")
# save_data("http://www.cherryblossom.org/results/2001/oof_f.html")

#save_data("http://www.cherryblossom.org/results/2002/ooff.htm")
#save_data("http://www.cherryblossom.org/results/2002/oofm.htm")

#save_data("http://www.cherryblossom.org/results/2003/CB03-F.HTM")
#save_data("http://www.cherryblossom.org/results/2003/CB03-M.HTM")

#save_data("http://www.cherryblossom.org/results/2004/women.htm")
#save_data("http://www.cherryblossom.org/results/2004/men.htm")

#save_data("http://www.cherryblossom.org/results/2005/CB05-F.htm")
#save_data("http://www.cherryblossom.org/results/2005/CB05-M.htm")

#save_data("http://www.cherryblossom.org/results/2006/women.htm")
#save_data("http://www.cherryblossom.org/results/2006/men.htm")

#save_data("http://www.cherryblossom.org/results/2007/women.htm")
#save_data("http://www.cherryblossom.org/results/2007/men.htm")

#save_data("http://www.cherryblossom.org/results/2008/women.htm")
#save_data("http://www.cherryblossom.org/results/2008/men.htm")

#save_data("http://www.cherryblossom.org/results/2009/09cucb-F.htm")
#save_data("http://www.cherryblossom.org/results/2009/09cucb-M.htm")

#save_data("http://www.cherryblossom.org/results/2010/2010cucb10m-f.htm")
#save_data("http://www.cherryblossom.org/results/2010/2010cucb10m-m.htm")

#save_data("http://www.cherryblossom.org/results/2011/2011cucb10m-f.htm")
#save_data("http://www.cherryblossom.org/results/2011/2011cucb10m-m.htm")

#save_data("http://www.cherryblossom.org/results/2012/2012cucb10m-f.htm")
#save_data("http://www.cherryblossom.org/results/2012/2012cucb10m-m.htm")

#save_data("http://www.cherryblossom.org/results/2012/2012cucb10m-f.htm")
#save_data("http://www.cherryblossom.org/results/2012/2012cucb10m-m.htm")

