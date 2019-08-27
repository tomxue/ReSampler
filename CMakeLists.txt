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

    if (WIN32)

        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            message(STATUS "Windows 32-bit")
            include_directories(libsndfile/include fftw32)
            link_directories(libsndfile32/lib fftw32)
        else()
            message(STATUS "Windows 64-bit")
            include_directories(libsndfile/include fftw64)
            link_directories(libsndfile/lib fftw64)
        endif()

        add_library(ReSamplerLib ${SOURCE_FILES})

        if(MSVC)
            message(STATUS "Visual Studio")
            target_link_libraries(ReSamplerLib libsndfile-1.lib libfftw3-3.lib ${CMAKE_THREAD_LIBS_INIT})

            if(NOT CMAKE_SIZEOF_VOID_P EQUAL 4)
                # create additional 64-bit config: ReleaseAVX
                set (CMAKE_CONFIGURATION_TYPES "Debug;Release;ReleaseAVX" CACHE STRING "" FORCE)

                # add /arch:AVX and /DUSE_AVX flags for ReleaseAVX configuration
                set (CMAKE_CXX_FLAGS_RELEASEAVX ${CMAKE_CXX_FLAGS_RELEASE_INIT} " /arch:AVX /DUSE_AVX")

                # same linker flags
                set (CMAKE_EXE_LINKER_FLAGS_RELEASEAVX ${CMAKE_EXE_LINKER_FLAGS_RELEASE_INIT})
            endif()

        else()
            message(STATUS "NOT Visual Studio")
            target_link_libraries(ReSamplerLib libsndfile-1 libfftw3-3 ${CMAKE_THREAD_LIBS_INIT})
        endif()

    else()

        if(APPLE)
            # since Mojave, cmake no longer finds these automatically:
            # (see https://stackoverflow.com/questions/54068035/linking-not-working-in-homebrews-cmake-since-mojave)
            include_directories(/usr/local/include)
            link_directories(/usr/local/lib)
        endif()

        add_library(ReSamplerLib SHARED ${SOURCE_FILES})
        target_link_libraries(ReSamplerLib fftw3 sndfile ${CMAKE_THREAD_LIBS_INIT})

    endif()

    add_executable(ReSampler main.cpp)
    target_link_libraries(ReSampler ReSamplerLib)

endif()
