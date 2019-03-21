#!/bin/sh
CFLDIR=$HOME/TurbPlasma/CFLIB
LIBDIR=$CFLDIR/lib
F90DIR=$CFLDIR/F90
if [ ! $1 ]; then
echo "Which Machine?"
read MACH
else 
export MACH=$1
fi

cd $F90DIR

if [[ "$MACH" =~ ^(flare) ]]; then
f2py3 $(\ls *.f90) -lfftw3 -c -m faf 
fi

if [[ "$MACH" =~ ^(goose) ]]; then
FFTWDIR=/home/tulasi/share/FFTW3/lib
f2py3 $(\ls *.f90) -L$FFTWDIR -lfftw3 -c -m faf 
fi

if [[ "$MACH" =~ ^(darter) ]]; then
module swap PrgEnv-intel PrgEnv-gnu
FFTWDIR=/opt/fftw/3.3.4.0/x86_64/lib
f2py3 $(\ls *.f90) -L$FFTWDIR -lfftw3 -c -m faf
fi

if [[ "$MACH" =~ ^(yellowstone) ]]; then
module load intel
FFTWDIR=/glade/p/work/tulasi/local/lib
f2py3 --opt=-free --fcompiler=intelem --f90flags="-O3" --compiler=intelem -L$FFTWDIR -lfftw3 -c -m faf $(\ls *.f90)
#f2py3 --debug-capi --opt=-free --fcompiler=intelem --compiler=intelem --f90flags="-check all -g -debug all -warn all -stand f08 -traceback" -L$FFTWDIR -lfftw3 -c -m faf $(\ls *.f90)
rm -f *__genmod.*
fi

if [[ "$MACH" =~ ^(cheyenne) ]]; then
ml mpt
ml fftw
f2py3 --build-dir $SCRATCH/f2py3/ --opt=-free --fcompiler=intelem --f90flags="-O3" --compiler=intelem -L/glade/p/work/tulasi/local/lib -lfftw3 -c -m faf $(\ls *.f90)
#FFTWDIR=/glade/p/work/tulasi/local/lib
#f2py3 --opt=-free --fcompiler=intelem --f90flags="-O3" --compiler=intelem -L$FFTWDIR -lfftw3 -c -m fafC $(\ls *.f90)
#f2py3 --debug-capi --opt=-free --fcompiler=intelem --compiler=intelem --f90flags="-check all -g -debug all -warn all -stand f08 -traceback" -L$FFTWDIR -lfftw3 -c -m faf $(\ls *.f90)
rm -f *__genmod.*
fi

##echo "f2py3 $FFTWDIR -lfftw3 -c -m faf $(\ls *.f90)"
#mv faf*.so ../lib/.
