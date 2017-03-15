/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef OSSPECIFIC_H
#define OSSPECIFIC_H 1

// macros which address differences between operating systems go here.

#define USE_CHRONO // use the std <chrono> library for timing.  

#ifdef _WIN32
#define NOMINMAX // disable min() and max() macros (use std:: library instead)

#pragma warning(disable : 4996) // suppress pointless MS "deprecation" warnings
#pragma warning(disable : 4244) // suppress double-to-float warnings

// Timer macros:
#ifndef USE_CHRONO
// don't use the C++11 <chrono> api ...
// Background: Older versions of MSVC prior to 2015 had issues with low timer resolution.

#include <windows.h>
#define START_TIMER() LARGE_INTEGER starttime,finishtime,elapsed,frequency,timetaken; \
	QueryPerformanceFrequency(&frequency); \
	QueryPerformanceCounter(&starttime)

#define STOP_TIMER() QueryPerformanceCounter(&finishtime); \
	elapsed.QuadPart=finishtime.QuadPart-starttime.QuadPart; \
	timetaken.QuadPart=((1000*elapsed.QuadPart)/frequency.QuadPart); \
	std::cout << "Time=" << static_cast<long>(timetaken.QuadPart) << " ms" << std::endl
#else
#include <chrono>
#define START_TIMER() auto beginTimer = std::chrono::high_resolution_clock::now()
#define STOP_TIMER() auto endTimer = std::chrono::high_resolution_clock::now(); \
std::cout << "Time=" << std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - beginTimer).count() << " ms" << std::endl
#endif // !USE_CHRONO

#else // Non-Windows:
typedef uint64_t __int64;
#define stricmp strcasecmp
// Timer macros:
#include <chrono>
#define START_TIMER() auto beginTimer = std::chrono::high_resolution_clock::now()
#define STOP_TIMER() auto endTimer = std::chrono::high_resolution_clock::now(); \
std::cout << "Time=" << std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - beginTimer).count() << " ms" << std::endl
#endif // ends Non-Windows

#endif // !OSSPECIFIC_H