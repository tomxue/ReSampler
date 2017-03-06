#ifndef OSSPECIFIC_H
#define OSSPECIFIC_H 1

// macros which address differences between operating systems go here

#ifdef _WIN32
#include <windows.h>
#define ZERO_64 0i64
// Timer macros:
#define START_TIMER() LARGE_INTEGER starttime,finishtime,elapsed,frequency,timetaken; \
	QueryPerformanceFrequency(&frequency); \
	QueryPerformanceCounter(&starttime)

#define STOP_TIMER() QueryPerformanceCounter(&finishtime); \
	elapsed.QuadPart=finishtime.QuadPart-starttime.QuadPart; \
	timetaken.QuadPart=((1000*elapsed.QuadPart)/frequency.QuadPart); \
	std::cout << "Time=" << static_cast<long>(timetaken.QuadPart) << " ms" << std::endl

#else // Non-Windows:

typedef uint64_t __int64;
#define ZERO_64 0LL
#define stricmp strcasecmp
// Timer macros:
#define START_TIMER() std::cout << "(start timer)" << std::endl
#define STOP_TIMER()  std::cout << "(end timer)" << std::endl

#endif // ends Non-Windows

#endif // !OSSPECIFIC_H