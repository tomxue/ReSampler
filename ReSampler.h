#ifndef RESAMPLER_H
#define RESAMPLER_H 1

#include <Windows.h>

const std::string strUsage("usage: resampler.exe -i <inputfile> [-o <outputfile>] -r <samplerate>\n");

#define BUFFERSIZE 4096
#define FILTERSIZE_HUGE 32767
#define FILTERSIZE_MEDIUM 511

#define MAX_8BIT 127
#define MAX_16BIT 32767
#define MAX_16BIT_MINUS_10TH_OF_DB 32392
#define MAX_24BIT 8388607

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
bool Convert(const std::string & InputFilename, const std::string & OutputFilename, unsigned int OutputSampleRate, bool bNormalize);
template<typename FloatType> bool makeBigAssLPF(FloatType* filter, int windowLength, FloatType transFreq, FloatType sampFreq);

// Timer macros:
#define START_TIMER() LARGE_INTEGER starttime,finishtime,elapsed,frequency,timetaken; \
	QueryPerformanceFrequency(&frequency); \
	QueryPerformanceCounter(&starttime)

#define STOP_TIMER() QueryPerformanceCounter(&finishtime); \
	elapsed.QuadPart=finishtime.QuadPart-starttime.QuadPart; \
	timetaken.QuadPart=((1000*elapsed.QuadPart)/frequency.QuadPart); \
	std::cout << "Time=" << static_cast<long>(timetaken.QuadPart) << " ms" << std::endl

#endif // !RESAMPLER_H
