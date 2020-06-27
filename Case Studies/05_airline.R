
x <- read.csv("1987.csv")
nrow(x)



totalSat <- 0
for (year in 1987:1988) {
   x <- read.csv(paste(year, ".csv", sep=''))
   totalSat <- totalSat + sum(x$DayOfWeek == 6)
}
totalSat



counts <- sapply(sprintf("%d.csv", 1987:1988),
                 function(f)
                    sum(read.csv(f)$DayOfWeek == 6))
sum(counts)

totalSat <- 0
for (year in 1987:1988) {
      x <- read.csv(paste(year, '.csv', sep=''))
      totalSat <- totalSat + sum(x$DayOfWeek == 6)
      rm(x)
      gc()
}
totalSat



library(RSQLite)
delay.con <- dbConnect("SQLite", dbname = "AirlineDelay.sqlite3")

delays87 <- dbGetQuery(delay.con, 
                   "SELECT * FROM AirlineDelay WHERE Year=1987")

nrow(delays87)


dbGetQuery(delay.con, "SELECT COUNT(*), Year FROM AirlineDelay
                                              WHERE Year=1987")


dbGetQuery(delay.con, 
          "SELECT COUNT(*), Year FROM AirlineDelay GROUP BY Year")


dbGetQuery(delay.con, 
           "SELECT COUNT(*), DayOfWeek FROM AirlineDelay
                                       WHERE DayOfWeek = 6")


x <- read.big.matrix("airline.csv", header = TRUE, 
                     backingfile = "airline.bin",
                     descriptorfile = "airline.desc",
                     type = "integer", extraCols = "age")

dim(x)  # How big is x?

x[1:6,1:6] # Show the first 6 rows and columns.


sum(x[, "Year"] == 1987)


sum(x[,"DayOfWeek"] == 6)


y <- attach.big.matrix("airline.desc")

foo <- big.matrix(nrow = 3, ncol = 3, type = "integer", init = 0)

foo


bar <- foo

bar[1,1] <- 1
foo


x <- attach.big.matrix("airline.desc")
dayCount = integer(7)
for (i in 1:7) 
  dayCount[i] <-  sum(x[,"DayOfWeek"] == i)

dayCount



state <- numeric(10)
for (i in 2:10) 
  state[i] <- state[i - 1] + sample( c(-1, 1), 1 )
state



library(foreach)
dayCount <- foreach(i = 1:7, .combine=c) %do% {
                                         sum(x[,"DayOfWeek"] == i)
                                       }

     # Split the rows of x by days of the week.
dow <- split(1:nrow(x), x[,"DayOfWeek"])
     # Rename the names of dow
names(dow) <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
     # Get the first 6 rows corresponding to Monday flights.
dow$Mon[1:6]


dayCount <- foreach(dayInds = dow, .combine = c) %do% {
                                                  length(dayInds)
                                                 }
dayCount



     # Divide CRSDepTime by 100 and take the floor to 
     # get the departure hour.
depHours <- floor(x[,"CRSDepTime"]/100)
     # Set the departure hours listed as 24 to 0.
depHours[depHours==24] <- 0

    # Split on the hours.
hourInds <- split(1:length(depHours), depHours)

    # Create a variable to hold the quantile probabilities.
myProbs <- c(0.9, 0.99, 0.999, 0.9999)

    # Use foreach to find the quantiles for each hour.
delayQuantiles <- foreach( hour = hourInds, .combine=cbind) %do% {
                     require(bigmemory)
                     x <- attach.big.matrix("airline.desc")
                     quantile(x[hour, "DepDelay"], myProbs, 
                                                   na.rm = TRUE)
                  }

    # Clean up the column names.
colnames(delayQuantiles) <- names(hourInds)

    # Load the parallel package so we can find 
    # how many cores are on the machine.
library(parallel)

    # Load our parallel backend.
library(doSNOW)

    # Use the total number of cores on the 
    # machine minus one.
numParallelCores <- max(1, detectCores()-1)

    # Create the parallel processes.
cl <- makeCluster(rep("localhost", numParallelCores), 
                     type = "SOCK")

    # Register the parallel processes with foreach.
registerDoSNOW(cl)

    # Run the foreach loop again, this time 
    # with %dopar% so that it is executed in parallel.
delayQuantiles <- foreach(hour=hourInds, .combine=cbind) %dopar% {
                require(bigmemory)
                x <- attach.big.matrix("airline.desc")
                quantile(x[hour, "DepDelay"], myProbs, na.rm=TRUE)
              }
colnames(delayQuantiles) <- names(hourInds)
stopCluster(cl)

library(ggplot2)
dq <- melt(delayQuantiles)
names(dq) <- c("percentile", "hour", "delay")
qplot(hour, delay, data = dq, color = percentile, geom = "line")

length(unique(x[,"TailNum"]))


planeStart <- foreach(tailInds = tailSplit, .combine=c) %dopar% {
        require(bigmemory)
        x <- attach.big.matrix("airline.desc")

            # Get the first year this tail code appears in the
            # data set.
        minYear <- min(x[tailInds, "Year"], na.rm = TRUE)
 
            # Get the rows that have the same year.
        minYearTailInds <- 
                tailInds[which(x[tailInds, "Year"] == minYear)]

            # The first month this tail code appears is the 
            # minimum month for rows indexed by minYearTailInds.
        minMonth <- min(x[minYearTailInds, "Month"], na.rm = TRUE)

            # Return the first time the tail code appears
            # in months A.D.
        12*minYear + minMonth
    }

x[,"age"] <- x[,"Year"] * 12 + x[,"Month"] -
                                    planeStart[x[,"TailNum"]]

library(biganalytics)
blm <- biglm.big.matrix( ArrDelay ~ age, data = x )
