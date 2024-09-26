# TopMGQuant
(global version)

## System requirements
GCC version higher than 4.8.2 for C++11 support

CMake (>= 3.1)

## Linux(Ubuntu)
    # install compiling tools
    sudo apt-get install build-essential cmake

    # install other library packages
    sudo apt-get install zlib1g-dev libboost-filesystem-dev \
                     libboost-program-options-dev \
                     libboost-system-dev \
                     libboost-thread-dev \
                     libboost-iostreams-dev \
                     libboost-chrono-dev \
                     libxalan-c-dev

    # install the catch unit test framework (https://github.com/philsquared/Catch)
    sudo apt-get install catch

    # install Qt5 for GUI
    sudo apt-get install qtbase5-dev

    # cd the toppic suite source folder toppic-suite-1.x.x
    # replace 1.x.x with the version number
    cd toppic-suite-1.x.x

    # build 
    mkdir build
    cd build
    cmake ..
    make -j$(nproc)

    # add the folder toppic_resources to the folder toppic_suite_1.x.x/bin
    cd ../bin
    ln -s ../toppic_resources .


    # test data
    protein: uniprot-st.fasta
    spectrum: st_1_ms1.feature  st_1_ms2.feature  st_1_ms2.msalign
    modification: variable_mods.txt

