#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/lib/x86_64-linux-gnu -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/opt/OpenBLAS/include -I/usr/include/x86_64-linux-gnu
OPTCLOCAL=-fPIC -march=native