# Compiling SIMOLANT from sources

## Linux

You need gcc/g++ and FLTK (V1.3) installed.

* For deb-based systems, use:<br />
  `sudo apt install g++ libfltk1.3-dev`
* Then, compile by:<br />
  `./make.sh`<br />
  or<br />
  `fltk-config --use-images --compile simolant.cc`

## Microsoft Windows

You need MinGW (with gcc/g++) and FLTK (V1.3) installed.

Warning: may be obsolete, I have one ancient installation working and could not get a new one on Windows 11

* For a dynamically linked version, see linux above.
* For statically-linked versions, both standard (windows) and console, use:<br />
  `sh winmake.sh`
* Hints:
  Windows desktop as seen from MinGW console looks like:<br />
  `/c/Windows/Users/<USER>/Desktop`<br />
  MinGW home seen from Windows may be:<br />
  `C:\MinGW\msys\1.0\home\<USER>`

## macOS

I am not a Mac user. I suspect that it is in not possible to run a home-made executable on a different machine without being approved by authorities

* To compile under macOS, Clang is recommended. You need FLTK (V1.3).
* Info on [the old page](http://old.vscht.cz/fch/software/simolant/index-en.html) may help
