CURRENT_DIR=`pwd`
INCLUDE_DIR=$CURRENT_DIR/../include
echo $INCLUDE_DIR
gcc -O3 -mfpmath=sse -finline-limit=5000 -fno-strict-aliasing  -msse2 -lboost_filesystem -lboost_program_options -lboost_regex -lboost_thread -I$INCLUDE_DIR __moldefs.cpp bmmkern.cpp -o bmmkern.exe
