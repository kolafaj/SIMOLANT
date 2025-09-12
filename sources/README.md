# Compiling SIMOLANT from sources

## Linux

You need gcc/g++ and FLTK (V1.1+).

* To install the compiler and FLTK (V1.3) in deb-based distros:<br />
  `sudo apt install g++ libx11-dev libfltk1.3-dev`
* Unzip SIMOLANT and select the `sources/` directory; e.g.:<br />
  `cd simolant/souces/`
* Compile SIMOLANT (with dynamic linking, optimized, debug flag, safe stack etc.):<br />
  `fltk-config --use-images --compile simolant.cc`

See also `make.sh`.

## Microsoft Windows

You need MinGW (linux-like terminal with gcc/g++) and FLTK (V1.1+).

### Install MinGW and FLTK

This guide was provided by ChatGPT and tested on Windows 11.

* Download and install MSYS2 (default settings are fine) from [msys2.org](https://www.msys2.org).
* Open the MSYS2 terminal (new blue icon with M2) and update:<br />
  `pacman -Syu`
* Use the 64-bit MinGW environment MSYS2 MINGW64 and install from the command line:<br />
  `pacman -S --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-fltk`<br />
  `pacman -S mingw-w64-x86_64-fltk mingw-w64-x86_64-fltk-images`
* <font color=gray>Optionally, use the (recommended) environment MSYS2 UCRT64:<br />
  `pacman -S --needed mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-fltk`<br />
  `pacman -S mingw-w64-ucrt-x86_64-fltk mingw-w64-ucrt-x86_64-fltk-images`</font>
* Verify:  <br />
  `g++ -dumpmachine' # should print: x86_64-w64-mingw32`<br />
  `fltk-config --version`

### Notes about the MSYS2/Windows interface:

* From the MSYS2 terminal, the Windows drives C:\\ D:\\ are visible as /c/ /d/.
* From the MSYS2 terminal, the Windows desktop is something like this (changes based on localization etc. are possible):<br />
  `/c/Windows/Users/<USER>/Desktop/`
* From the Windows GUI, the MSYS2 is something like this (depending on your setup):<br />
  `This PC → Windows (C:) → msys64 → home → <USER>`

### Compile SIMOLANT:

* Unzip SIMOLANT and select the `sources/` directory; e.g.:<br />
  `cd simolant/souces/`
* Compile SIMOLANT (simple, dynamic linking, to start from MSYS2):<br />
  `fltk-config --use-images --compile simolant.cc`
* As above with better control:<br />
  `LINK=-lfltk-images`<br />
  `g++ -O2 simolant.cc -o simolant.exe $(fltk-config --cxxflags --use-images) $LINK`
* Compile with partly static linking:<br />
  `LINK="-static-libgcc -static-libstdc++ -static -lfltk-images"`<br />
  `g++ -O2 simolant.cc -o simolant.exe $(fltk-config --ldstaticflags --use-images) $LINK`<br />
  The following libraries are not static and must be added to the installation ZIP:<br />
  `cd /c/msys64/mingw64/bin/` # or similar<br />
  `cp zlib1.dll libwinpthread-1.dll libpng16-16.dll libjpeg-8.dll libgcc_s_seh-1.dll DEST/`<br />
  where DEST/ is the directory with simolant.exe and simolant.html

See also `winmake.sh`.

## macOS

I am not a Mac user. I suspect that it is in not possible to run a home-made
executable on a different machine without being approved by authorities.  Bases on info from my former students:

* To compile under macOS, Clang is recommended. You need FLTK (V1.1+).
* Info on [the old page](http://old.vscht.cz/fch/software/simolant/index-en.html) may help
