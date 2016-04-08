#!/bin/bash
#$ -V
#$ -cwd
#$ -pe mpi 64

DELTA=0.0000001
ITER=2500

make clean
make parallel
rm -rf out
mkdir -p out

for t in 1 2 4 8 16 32 64; do
	mpiexec -n $t ./reconstruct.parallel -i $ITER -d $DELTA images/edgenew1024x1280.pgm -o out/1024x1280-$t.pgm
	mpiexec -n $t ./reconstruct.parallel -i $ITER -d $DELTA images/edgenew512x384.pgm -o out/512x384-$t.pgm
	mpiexec -n $t ./reconstruct.parallel -i $ITER -d $DELTA images/edgenew768x768.pgm -o out/768x768-$t.pgm
	mpiexec -n $t ./reconstruct.parallel -i $ITER -d $DELTA images/edgenew192x128.pgm -o out/192x128-$t.pgm
	mpiexec -n $t ./reconstruct.parallel -i $ITER -d $DELTA images/edgenew256x192.pgm -o out/256x192-$t.pgm
done

make clean
make serial
./reconstruct.serial -i $ITER -d $DELTA images/edgenew1024x1280.pgm -o out/1024x1280-serial.pgm
./reconstruct.serial -i $ITER -d $DELTA images/edgenew512x384.pgm -o out/512x384-serial.pgm
./reconstruct.serial -i $ITER -d $DELTA images/edgenew768x768.pgm -o out/768x768-serial.pgm
./reconstruct.serial -i $ITER -d $DELTA images/edgenew192x128.pgm -o out/192x128-serial.pgm
./reconstruct.serial -i $ITER -d $DELTA images/edgenew256x192.pgm -o out/256x192-serial.pgm

sha256sum -c sha256-checksum
