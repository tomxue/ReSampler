cmake_minimum_required(VERSION 3.6)
project(ReSampler)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -pthread -m64")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(WIN32)

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
            fraction.h
            noiseshape.h
            osspecific.h
            raiitimer.h
            ReSampler.cpp
            ReSampler.h
            srconvert.h)

    include_directories(libsndfile/include fftw64)
    link_directories(libsndfile/lib fftw64)
    link_libraries(libsndfile-1  libfftw3-3)
    add_executable(ReSampler ${SOURCE_FILES})

elseif (APPLE)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
    link_libraries(-lfftw3 -lsndfile)

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
            fraction.h
            noiseshape.h
            osspecific.h
            raiitimer.h
            ReSampler.cpp
            ReSampler.h
            srconvert.h)

    add_executable(ReSampler ${SOURCE_FILES})
elseif(ANDROID)

    set(SOURCE_FILES
            libsndfilearm64_v8a/debug/src/sndfile.h
            libsndfilearm64_v8a/debug/src/sndfile.hh
            fftwarm64_v8a/debug/fftw3.h
            alignedmalloc.h
            biquad.h
            conversioninfo.h
            dff.h
            ditherer.h
            dsf.h
            FIRFilter.h
            fraction.h
            noiseshape.h
            osspecific.h
            raiitimer.h
            ReSampler.cpp
            ReSampler.h
            srconvert.h)

    include_directories(libsndfilearm64_v8a/debug/src fftwarm64_v8a/debug)
    link_directories(libsndfilearm64_v8a/debug fftwarm64_v8a/debug)
    add_library(ReSampler SHARED ${SOURCE_FILES})
    target_link_libraries(ReSampler sndfile fftw3)
else ()
    include_directories(/usr/include)
    link_directories(/usr/lib)
    link_libraries(-lfftw3 -lsndfile)

    set(SOURCE_FILES
            /usr/include/sndfile.h
            /usr/include/sndfile.hh
            /usr/include/fftw3.h
            alignedmalloc.h
            biquad.h
            conversioninfo.h
            dff.h
            ditherer.h
            dsf.h
            FIRFilter.h
            fraction.h
            noiseshape.h
            osspecific.h
            raiitimer.h
            ReSampler.cpp
            ReSampler.h
            srconvert.h csv.h)

    add_executable(ReSampler ${SOURCE_FILES})
endif()
