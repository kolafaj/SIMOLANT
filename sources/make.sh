#!/bin/bash
if [ $# -lt 1 ]; then
  cat <<EOF
Compile simolant. Call by:
  make.sh {s|d|D}
where
  g = -g (debug)
  s = static linking
  d = dynamic linking
  D = dynamic linking, move to ~/bin/ on success
EOF
  exit 1
fi

export LC_ALL=C

case $1 in
  g ) g++ -I/usr/include/freetype2 -g -D_THREAD_SAFE -D_REENTRANT -Wall -o simolant simolant.cc -lfltk_images -lfltk ;;
  d ) g++ -I/usr/include/freetype2 -O2 -D_THREAD_SAFE -D_REENTRANT -Wall -o simolant simolant.cc -lfltk_images -lfltk ;;
  D ) g++ -I/usr/include/freetype2 -O2 -D_THREAD_SAFE -D_REENTRANT -Wall -o simolant simolant.cc -lfltk_images -lfltk && mv simolant ~/bin/ ;;
  s ) fltk-config --use-images --compile simolant.cc ;;
  * ) echo "wrong argument, use one of: d D s"
esac
