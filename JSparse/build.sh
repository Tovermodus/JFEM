javac -h src/native -d out/java src/com/simonschmidt/*
echo "JAVA COMPILED"
g++ -llapack -lopenblas -c -fPIC -I/usr/lib/jvm/java-11-openjdk/include -I/usr/lib/jvm/java-11-openjdk/include/linux src/native/sparseSolver.cpp -o out/native/sparseSolver.o
echo "C++ COMPILED"
g++ -llapack -lopenblas -shared -fPIC -o out/native/libnative.so out/native/sparseSolver.o -lc
echo "BUILT SO"
