mkdir 3rd_party
cd 3rd_party

# KMC
git clone https://github.com/refresh-bio/KMC.git
cd KMC
make
cd ../


# STXXL
git clone https://github.com/stxxl/stxxl.git
cd stxxl
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install -DBUILD_STATIC_LIBS=ON
make
make install
cd ../../


# SDSL
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite/
/usr/bin/time sh install.sh install
cd ../


# SPARSEPP
git clone https://github.com/greg7mdp/sparsepp.git

# Google Test
git clone https://github.com/google/googletest.git
cd googletest/googletest
mkdir build
cd build 
cmake ../
make

