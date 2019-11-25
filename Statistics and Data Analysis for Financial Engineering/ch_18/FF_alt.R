# Stock/Bond/FX data.
stocks <- as.data.table(read.csv(paste0(data.dir, "Stock_FX_Bond_2004_to_2005.csv"), 
                                 header=T))
stocks$Date <- as.Date(stocks$Date, format = "%d-%b-%y")
stocks_subset <- stocks[, .(Date, GM_AC, F_AC, UTX_AC, MRK_AC)]
stocks_diff <- data.table(Date = stocks_subset[-1]$Date, 
                          apply(log(stocks_subset[, .(GM_AC, F_AC, UTX_AC, MRK_AC)]), 2, diff))

# Fama-French data.
FF_data <- as.data.table(read.table(paste0(data.dir, "FamaFrenchDaily.txt"), 
                                    header=T))

stocks = read.csv(paste0(data.dir, "Stock_FX_Bond_2004_to_2005.csv"), header=T)

FF_data = read.table(paste0(data.dir, "FamaFrenchDaily.txt"), header=T)

start_index = which( FF_data$date == "20040102" )
end_index = which( FF_data$date == "20051230" )
FF_data = FF_data[ start_index:end_index, ]
FF_data = FF_data[-1,]     # delete first row since stocks_diff loses a row due to differencing

start_index = which( stocks$DATE == "1/2/2004" )
end_index = which( stocks$DATE == "12/30/2005" )
stocks = stocks[ start_index:end_index, ]

stocks_subset = as.data.frame(cbind(stocks$GM_AC,stocks$F_AC,stocks$UTX_AC,stocks$MRK_AC))

stocks_diff = as.data.frame(100*apply(log(stocks_subset), 2, diff) - FF_data$RF) # Excess returns

names(stocks_diff) = c("GM", "Ford", "UTX", "Merck")

fit2=lm(as.matrix(stocks_diff)~FF_data$Mkt.RF + FF_data$SMB + FF_data$HML)
options(digits = 4)

coef(fit2)
