#ifndef RESAMPLER_H
#define RESAMPLER_H 1

#include <Windows.h>

const std::string strUsage("usage: resampler.exe -i <inputfile> [-o <outputfile>] -r <samplerate>\n");

#define BUFFERSIZE 8192


typedef struct fraction {
	int numerator;
	int denominator;
} Fraction;

int gcd(int a, int b);
Fraction GetSimplifiedFraction(int InputSampleRate, int OutputSampleRate);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, std::string & Parameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, unsigned int & nParameter);
bool findCmdlineOption(char ** begin, char ** end, const std::string & option);
template<typename FloatType>
bool Convert(const std::string & InputFilename, const std::string & OutputFilename, unsigned int OutputSampleRate, FloatType Limit);
template<typename FloatType> bool makeLPF(FloatType* filter, int Length, FloatType transFreq, FloatType sampFreq);

template<typename FloatType>
bool applyBlackmanWindow(FloatType * filter, int Length);

template<typename FloatType>
bool applyKaiserWindow(FloatType * filter, int Length, FloatType Beta);

template<typename FloatType>
FloatType I0(FloatType z); // 0th-order Modified Bessel function of the first kind

// Timer macros:
#define START_TIMER() LARGE_INTEGER starttime,finishtime,elapsed,frequency,timetaken; \
	QueryPerformanceFrequency(&frequency); \
	QueryPerformanceCounter(&starttime)

#define STOP_TIMER() QueryPerformanceCounter(&finishtime); \
	elapsed.QuadPart=finishtime.QuadPart-starttime.QuadPart; \
	timetaken.QuadPart=((1000*elapsed.QuadPart)/frequency.QuadPart); \
	std::cout << "Time=" << static_cast<long>(timetaken.QuadPart) << " ms" << std::endl

#endif // !RESAMPLER_H
