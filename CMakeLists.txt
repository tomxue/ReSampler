cmake_minimum_required(VERSION 3.8)
project(ReSampler)

#g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -O3 -L/usr/local/lib -I/usr/local/include

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
link_libraries(-lfftw3 -lsndfile)

add_compile_options(-pthread -O3)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        /usr/local/include/sndfile.h
        /usr/local/include/sndfile.hh
        /usr/local/include/fftw3.h
        alignedmalloc.h
        biquad.h
        conversioninfo.h
        dff.h
        ditherer.h
        dsf.h
        FIRFilter.h
        FIRFilterAVX.h
        fraction.h
        noiseshape.h
        osspecific.h
        raiitimer.h
        ReSampler.cpp
        ReSampler.h
        srconvert.h)

add_executable(ReSampler ${SOURCE_FILES})