#!/bin/sh

g++ -I/usr/include/freetype2 -O2 -D_THREAD_SAFE -D_REENTRANT -Wall -o simolant simolant.cc -lfltk_images -lfltk 

# make for Linux (FLTK), dynamic linking
# fltk-config --use-images --compile simolant.cc
