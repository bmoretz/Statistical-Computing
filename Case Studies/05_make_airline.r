library(iotools)

col_types <- c(rep("integer", 8), "character", "integer", "character",
  rep("integer", 5), "character", "character",
  rep("integer", 4), "character", rep("integer", 6))

col_names <- c("Year", "Month", "DayofMonth", "DayOfWeek", "DepTime", 
  "CRSDepTime", "ArrTime", "CRSArrTime", "UniqueCarrier", "FlightNum", 
  "TailNum", "ActualElapsedTime", "CRSElapsedTime", "AirTime", "ArrDelay", 
  "DepDelay", "Origin", "Dest", "Distance", "TaxiIn", "TaxiOut", "Cancelled", 
  "CancellationCode", "Diverted", "CarrierDelay", "WeatherDelay", "NASDelay", 
  "SecurityDelay", "LateAircraftDelay")

# A mutable closure to skip the header line.
make_airline_df_gen <- function() {
  first_run <- TRUE
  function(chunk) {
    x <- dstrsplit(chunk, col_types=col_types, sep=",", 
      skip=as.integer(first_run))
    colnames(x) <- col_names
    first_run <<- FALSE
    x
  }
}

#######
# One pass to get the unique character types.
#######

make_airline_df <- make_airline_df_gen()

us <- chunk.apply("airline.csv",
    # Get the unique values for each of the character columns.
    function(chunk) {
      x <- make_airline_df(chunk)
      list(carrier=unique(x$UniqueCarrier),
           tail_num=unique(x$TailNum),
           origin=unique(x$Origin),
           dest=unique(x$Dest),
           cancel_code=unique(x$CancellationCode))
    }, CH.MERGE=list)

########
# Make the maps between the character values and an integer.
########

unique_reduce <- function(x, y) unique(c(x, y))

carrier_names <- 
  sort(Reduce(unique_reduce, Map(function(x) x$carrier, us)))
carrier <- 1:length(carrier_names)
names(carrier) <- carrier_names

tail_num_names <- sort(Reduce(unique_reduce, Map(function(x) x$tail_num, us)))
tail_num_names <- tail_num_names[tail_num_names != "NA"]
tail_num <- 1:length(tail_num_names)
names(tail_num) <- tail_num_names

origin_names <- sort(Reduce(unique_reduce, Map(function(x) x$origin, us)))
origin <- 1:length(origin_names)
names(origin) <- origin_names

dest_names <- sort(Reduce(unique_reduce, Map(function(x) x$dest, us)))
dest <- 1:length(dest_names)
names(dest) <- dest_names

cancel_names <- sort(Reduce(unique_reduce, Map(function(x) x$cancel_code, us)))
cancel_names <- cancel_names[cancel_names != "" & cancel_names != "NA"] 
cancel_code <- 1:length(cancel_names)
names(cancel_code) <- cancel_names

save(carrier, tail_num, origin, dest, cancel_code, 
  file="airline_character_maps.RData", ascii=TRUE)

#######
# Now create the output file. 
#######

out_file <- file("airline_int_cols.csv", "wb")
writeBin(as.output(matrix(col_names, nrow=1), sep=","), out_file)

# Make sure to create a new make_airline_df mutable closure. first_run
# is FALSE in the old one.
make_airline_df <- make_airline_df_gen()

chunk.apply("airline.csv",
  function(chunk) {
    x <- make_airline_df(chunk)
    x$UniqueCarrier <- as.vector(as.integer(carrier[x$UniqueCarrier]))
    x$TailNum <- as.vector(as.integer(carrier[x$TailNum]))
    x$Origin <- as.vector(as.integer(carrier[x$Origin]))
    x$Dest <- as.vector(as.integer(carrier[x$Dest]))
    x$CancellationCode <- as.vector(as.integer(carrier[x$CancellationCode]))
    writeBin(as.output(x, sep=","), out_file)
    NULL
  })

close(out_file)

