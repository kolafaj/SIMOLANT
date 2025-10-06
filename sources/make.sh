#!/bin/sh
if [ $# -lt 1 ]; then
  cat <<EOF
Compile simolant (linux) with dynamic linking. Call by:
  make.sh {c,d,g,D}
where
  c = safest, optimized, debug flag (via fltk-config --compile)
  g = not optimized, debug flag (use gdb)
  o = optimized, no debug info
  i = install, as above + move to ~/bin/ on success
EOF
  exit 1
fi

export LC_ALL=C

# ? -I/usr/include/freetype2

case $1 in
  c ) fltk-config --use-images --compile simolant.cc ;;
  g ) g++ -g -Wall -o simolant simolant.cc -lfltk_images -lfltk ;;
  o ) g++ -O2 -Wall -o simolant simolant.cc -lfltk_images -lfltk ;;
  i ) g++ -O2 -Wall -o simolant simolant.cc -lfltk_images -lfltk && mv simolant ~/bin/ ;;
  * ) echo "Wrong argument, use:"; echo "  make.sh {c,g,o,i}"
esac
