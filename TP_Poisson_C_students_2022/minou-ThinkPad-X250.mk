#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/lib/x86_64-linux-gnu -llapack -llapacke -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native