g++ -I/usr/local/include -I/usr/local/include/FL/images -mwindows -DWIN32 \
  -o simolant32 simolant.cc \
  -mwindows --static /usr/local/lib/libfltk_images.a \
  /usr/local/lib/libfltk_png.a \
  /usr/local/lib/libfltk_z.a \
  /usr/local/lib/libfltk_jpeg.a \
  /usr/local/lib/libfltk.a \
  -lole32 -luuid -lcomctl32
