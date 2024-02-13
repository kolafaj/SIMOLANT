Compiling SIMOLANT from sources
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Linux
  You need gcc/g++ and FLTK (V1.3) installed.
  For deb-based systems, use:
    sudo apt install g++ libfltk1.3-dev
  Compile by:
    ./make.sh
  or
    fltk-config --use-images --compile simolant.cc

Microsoft Windows
  Yu need MinGW (with gcc/g++) + FLTK (V1.3) installed.
  For a dynamically linked version, see linux above.
  For statically-linked versions, both standard (windows) and console, use:
    sh winmake.sh
  Hints:
  Windows desktop as seen from MinGW console looks like:
    /c/Windows/Users/<USER>/Desktop
  MinGW home seen from Windows may be:
    C:\MinGW\msys\1.0\home\<USER>

macOS
  To compile under macOS, Clang is recommended. You need FLTK (V1.3).
