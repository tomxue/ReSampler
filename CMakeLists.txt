cmake_minimum_required(VERSION 3.8)
project(ReSampler)

include_directories(libsndfile/include fftw64)
link_directories(libsndfile/lib fftw64)
link_libraries(-llibsndfile-1 -llibfftw3-3)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        libsndfile/include/sndfile.h
        libsndfile/include/sndfile.hh
        fftw64/fftw3.h
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