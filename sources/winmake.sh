#!/bin/sh
if [ $# -lt 1 ]; then
  cat <<EOF
Compile simolant (windows/MinGW/MSYS2). Call by:
  make.sh {c,s,d}
where
  s = mostly static linking, optimized, no debug info
      libraries needed (will be packed with):
        zlib1.dll libwinpthread-1.dll libpng16-16.dll
        libjpeg-8.dll libgcc_s_seh-1.dll
  c = safest, optimized, debug flag, dynamic linking
  d = optimized, no debug info, dynamic linking
EOF
  exit 1
fi

export LC_ALL=C

case $1 in
  c ) fltk-config --use-images --compile simolant.cc ;;
  s ) LINK="-static-libgcc -static-libstdc++ -static -lfltk-images"
      g++ -O2 simolant.cc -o simolant.exe $(fltk-config --ldstaticflags --use-images) $LINK ;;
  d ) LINK=-lfltk-images
      g++ -O2 simolant.cc -o simolant.exe $(fltk-config --cxxflags --use-images) $LINK ;;
  * ) echo "Wrong argument, use:"; echo "  make.sh {c,s,d}"
esac

ls -l simolant.exe
