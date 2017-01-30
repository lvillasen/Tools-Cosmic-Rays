#################################################################################################################
#Modelo de JF2012, descrito en http://arxiv.org/abs/1204.3662
#################################################################################################################

# El archivo TA.dat.dat.GAL tiene este formato
# n theta Energy RA DEC GalLong GalLat  GalLong_B GalLat_B Delta

# El archivo input.txt tiene este formato
# Charge,Mass,Energy, galactic_l, galactic_b (in degrees)
# Usage sh Backprop.sh p
model='jf2012'
#mass='p'
#mass='O'
#mass='Fe'
mass=$1

file_size=`wc  Data.dat.GAL | awk '{print $1}'`

if [ -f "./t" ]; then rm ./t ; fi; 
for i in `eval echo {1..$file_size}`; 
do printf $i; printf " "  
if [ $mass = 'p' ]; then \
    awk -v LINE=${i} '(NR==LINE)\
    {printf "%s\t %s\t %5.3f\t %5.3f\t %5.3f\t \n","1","1", $3,$6,$7}' \
    ./Data.dat.GAL > ./input.txt; fi 
if [ $mass = 'O' ]; then \
    awk -v LINE=${i} '(NR==LINE)\
    {printf "%s\t %s\t %5.3f\t %5.3f\t %5.3f\t \n","8","16", $3,$6,$7}' \
    ./Data.dat.GAL > ./input.txt; fi 
if [ $mass = 'Fe' ]; then \
    awk -v LINE=${i} '(NR==LINE)\
    {printf "%s\t %s\t %5.3f\t %5.3f\t %5.3f\t \n","26","56", $3,$6,$7}' \
    ./Data.dat.GAL > ./input.txt; fi 
julia CRBackProp.jl | \
awk  '(NR==2){printf "%5.3f\t %5.3f %5.3f\n", $7,$8,$9}' >> ./t ;
done ; 
echo $model $mass
paste ./Data.dat.GAL ./t > ./Data.$mass.$model.dat
ls -alt Data*