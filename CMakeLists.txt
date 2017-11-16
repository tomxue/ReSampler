cmake_minimum_required(VERSION 3.8)
project(ReSampler)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
link_libraries(-lfftw3 -lsndfile)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto -DNDEBUG")

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