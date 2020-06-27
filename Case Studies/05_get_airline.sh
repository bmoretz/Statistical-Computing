# Get the airline data.
for year in {1987..2008}
do
  wget http://stat-computing.org/dataexpo/2009/$year.csv.bz2 
done

# Write the uncompressed data to a single file.
pbzip2 -dc 1987.csv.bz2 > airline.csv
for year in {1988..2008}
do
  pbzip2 -dc $year.csv.bz2 | sed -e "1d" >> airline.csv
done
