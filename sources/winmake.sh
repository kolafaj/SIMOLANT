
# make for Windows (MinGW+FLTK), static linking
CXXFLAGS=$(fltk-config --use-images --cxxflags)
LDSTATICFLAGS=$(fltk-config --use-images --ldstaticflags)

echo "CXXFLAGS:      $CXXFLAGS"
echo "LDSTATICFLAGS: $LDSTATICFLAGS"

# standard windows version (console, option -D does not work)
g++ ${CXXFLAGS} -o simolant simolant.cc ${LDSTATICFLAGS} -static-libgcc -static-libstdc++ -mwindows

# console version (debug option -D active if started from the MinGW console)
g++ ${CXXFLAGS} -o simolant-console simolant.cc ${LDSTATICFLAGS} -static-libgcc -static-libstdc++ -mconsole

