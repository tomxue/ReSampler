project(ReSampler C CXX)
cmake_minimum_required(VERSION 3.6)
find_package(Threads REQUIRED)
set(CMAKE_CXX_STANDARD 11)

if(ANDROID)

    add_subdirectory (android/src/fftw-3.3.8)
    set(BUILD_SHARED_LIBS_SAVED "${BUILD_SHARED_LIBS}")
    set(BUILD_PROGRAMS_SAVED "${BUILD_PROGRAMS}")
    set(BUILD_EXAMPLES_SAVED "${BUILD_EXAMPLES}")
    set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
    set(BUILD_SHARED_LIBS OFF)
    set(BUILD_PROGRAMS OFF)
    set(BUILD_EXAMPLES OFF)
    set(BUILD_TESTING OFF)
    add_subdirectory (android/src/libsndfile)
    set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_SAVED}")
    set(BUILD_PROGRAMS "${BUILD_PROGRAMS_SAVED}")
    set(BUILD_EXAMPLES "${BUILD_EXAMPLES_SAVED}")
    set(BUILD_TESTING "${BUILD_TESTING_SAVED}")
    include_directories(android/src/fftw-3.3.8/api)
    include_directories($CMAKE_SOURCE_DIR/android/src/libsndfile/src)

    set(SOURCE_FILES
        alignedmalloc.h
        biquad.h
        conversioninfo.h
        conversioninfo.cpp
        csv.h
        dff.h
        ditherer.h
        dsf.h
        FIRFilter.h
        fraction.h
        factorial.h
        noiseshape.h
        osspecific.h
        raiitimer.h
        main.cpp
        ReSampler.cpp
        ReSampler.h
        srconvert.h)

    add_library(ReSampler SHARED ${SOURCE_FILES})
    target_link_libraries(ReSampler sndfile fftw3 log)

elseif (WIN32)

    set(SOURCE_FILES
        alignedmalloc.h
        biquad.h
        conversioninfo.h
        csv.h
        dff.h
        ditherer.h
        dsf.h
        FIRFilter.h
        fraction.h
        factorial.h
        noiseshape.h
        osspecific.h
        raiitimer.h
        ReSampler.cpp
        ReSampler.h
        srconvert.h
        conversioninfo.cpp)

    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        message(STATUS "Windows 32-bit")
        include_directories(libsndfile/include fftw32)
        link_directories(libsndfile32/lib fftw32)
    else()
        message(STATUS "Windows 64-bit")
        include_directories(libsndfile/include fftw64)
        link_directories(libsndfile/lib fftw64)
    endif()


    if(MSVC_VERSION GREATER 1700)
        message(STATUS "Visual Studio")
        add_library(ReSamplerLib ${SOURCE_FILES})
        add_executable(ReSampler main.cpp)
        target_link_libraries(ReSamplerLib libsndfile-1.lib libfftw3-3.lib ${CMAKE_THREAD_LIBS_INIT})
        target_link_libraries(ReSampler ReSamplerLib)
    else()
        message(STATUS "NOT Visual Studio")
        add_library(ReSamplerLib SHARED ${SOURCE_FILES})
        target_link_libraries(ReSamplerLib libsndfile-1 libfftw3-3 ${CMAKE_THREAD_LIBS_INIT})
        add_executable(ReSampler main.cpp)
        target_link_libraries(ReSampler ReSamplerLib)
    endif()

else()
    set(SOURCE_FILES
        alignedmalloc.h
        biquad.h
        conversioninfo.h
        csv.h
        dff.h
        ditherer.h
        dsf.h
        FIRFilter.h
        fraction.h
        factorial.h
        noiseshape.h
        osspecific.h
        raiitimer.h
        ReSampler.cpp
        ReSampler.h
        srconvert.h
        conversioninfo.cpp)

    add_library(ReSamplerLib SHARED ${SOURCE_FILES})
    target_link_libraries(ReSamplerLib ReSamplerLib fftw3 sndfile ${CMAKE_THREAD_LIBS_INIT})
    add_executable(ReSampler main.cpp)
    target_link_libraries(ReSampler ReSamplerLib)

endif()
