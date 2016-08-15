#!/bin/bash 

gfortran constant.f90 nearest.f90 init.f90 solveit.f90 graphene.f90 -o graphene.exe -llapack

./graphene.exe


colm=2
echo 'input number of eigenenergies to view'
read num

max=$[2+$num] 

while [ $colm -lt $max ];

do A+=' -bxy 1:'$colm 
  let colm=colm+1
done

xmgrace -block eigen.dat  $A -param ar.par #&
#-param ar.par 
