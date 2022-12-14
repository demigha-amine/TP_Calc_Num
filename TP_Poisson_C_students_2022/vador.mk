#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/lib/x86_64-linux-gnu -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/include
OPTCLOCAL=-fPIC -march=native