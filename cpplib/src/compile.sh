#GPP="g++ -Wall -fbounds-check "
#GPP="g++ -I /home/jose/apps/eigen/"
GPP="g++ -O3 $(pkg-config --cflags eigen3)  -fbounds-check"
GPP="g++ -I/home/jose/apps/eigen-3.2.5  -fbounds-check -larmadillo"
GPP="g++ -fmem-report -Wall -I/home/jose/apps/eigen-3.2.5  -fbounds-check"
GPP="g++ -Wall -I/home/jose/apps/eigen-3.2.5  -fbounds-check"
#GPPF=$GPP" -c -g -eigen"

rm *.o *.gch
#$GPPF objects.h
#$GPPF functions.h
#$GPPF generate_basis.cpp
#$GPPF get_basis_nelectrons.cpp
#$GPPF main.cpp
$GPP -ggdb main.cpp -o main.x 

