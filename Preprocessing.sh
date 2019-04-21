Obs=$1
if [ "$Obs" = "TA" ]; then
echo "Processing data from the $Obs Observatory .."
awk '(substr($0,1,1) != "#"){n+=1;printf n " " $5 " " $6 " " $7 " " $8 "\n"} ' TA_AstrophysJLett_2014.dat \
|  sed 's/−/-/g' > Data.dat; fi
if  [ "$Obs" = "Both" ]; then
echo "Processing data from  $Obs Observatories .."
awk '(substr($0,1,1) != "#"){n+=1;printf n " " $5 " " $6 " " $7 " " $8 "\n"} ' TA_AstrophysJLett_2014.dat \
|  sed 's/−/-/g' > Data.dat ;
awk '(substr($0,1,1) != "#"){n+=1;printf n " " $3 " " $4 " " $5 " " $6 "\n"} ' Auger_AstrophysJ804_2015.dat \
|  sed 's/−/-/g' >> Data.dat; fi
if  [ "$Obs" = "Auger" ]; then
echo "Processing data from the $Obs Observatory .."
awk '(substr($0,1,1) != "#"){n+=1;printf n " " $3 " " $4 " " $5 " " $6 "\n"} ' Auger_AstrophysJ804_2015.dat \
|  sed 's/−/-/g' > Data.dat; fi

echo "n Theta(º) E(EeV) RA(º) DEC(º)"
head Data.dat
echo "Number of events";
wc  Data.dat | awk '{print $1}'
