javac -h src/native -d out/java src/com/simonschmidt/*
g++ -llapacl -lblas -c -fPIC -I/usr/lib/jvm/java-11-openjdk/include -I/usr/lib/jvm/java-11-openjdk/include/linux src/native/sparseSolver.cpp -o out/native/sparseSolver.o
g++ -llapack -lblas -shared -fPIC -o out/native/libnative.so out/native/sparseSolver.o -lc