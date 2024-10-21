#!/bin/bash
if [ $# -lt 1 ]; then
  cat <<EOF
Compile simolant. Call by:
  make.sh {s|d}
where
  s = static linking
  d = dynamic linking
EOF
  exit 1
fi

export LC_ALL=C

case ${1:0:1} in
  d | D ) g++ -I/usr/include/freetype2 -O2 -D_THREAD_SAFE -D_REENTRANT -Wall -o simolant simolant.cc -lfltk_images -lfltk ;;
  s | S ) fltk-config --use-images --compile simolant.cc ;;
  * ) echo "wrong argument"
esac
