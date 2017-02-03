#ifndef RESAMPLER_H
#define RESAMPLER_H 1

#include <map>
#include <Windows.h>
#include <sndfile.h>
#include <sndfile.hh>

const std::string strVersion("1.1.7-pre-release(NOT FOR PRODUCTION !)");
const std::string strUsage("usage: resampler.exe -i <inputfile> [-o <outputfile>] -r <samplerate> [-b <bitformat>] [-n [<normalization factor>]]\n");
const std::string strExtraOptions("--help\n--version\n--doubleprecision\n--listsubformats <ext>\n--dither [<amount>] [--autoblank]\n--minphase\n--flacCompression <compressionlevel>\n--vorbisQuality <quality>\n--noClippingProtection\n");

#define BUFFERSIZE 8192 // buffer size for file reads

#pragma warning(disable : 4996) // suppress pointless MS "deprecation" warnings
#pragma warning(disable : 4244) // suppress double-to-float warnings

typedef struct fraction {
	int numerator;
	int denominator;
} Fraction;

// map of commandline subformats to libsndfile subformats:
const std::map<std::string,int> subFormats = { 
	{ "s8",SF_FORMAT_PCM_S8 },
	{ "u8",SF_FORMAT_PCM_U8 },
	{ "8",SF_FORMAT_PCM_U8 },	// signed or unsigned depends on major format of output file eg. wav files unsigned
	{ "16", SF_FORMAT_PCM_16 },
	{ "24", SF_FORMAT_PCM_24 },
	{ "32", SF_FORMAT_PCM_32 },
	{ "32f",SF_FORMAT_FLOAT },
	{ "64f",SF_FORMAT_DOUBLE },
	{ "ulaw",SF_FORMAT_ULAW },
	{ "alaw",SF_FORMAT_ALAW },
	{ "ima-adpcm",SF_FORMAT_IMA_ADPCM },
	{ "ms-adpcm",SF_FORMAT_MS_ADPCM },
	{ "gsm610",SF_FORMAT_GSM610 },
	{ "vox-adpcm",SF_FORMAT_VOX_ADPCM },
	{ "g721-32",SF_FORMAT_G721_32 },
	{ "g723-24",SF_FORMAT_G723_24 },
	{ "g723-40",SF_FORMAT_G723_40 },
	{ "dwvw12",SF_FORMAT_DWVW_12 },
	{ "dwvw16",SF_FORMAT_DWVW_16 },
	{ "dwvw24",SF_FORMAT_DWVW_24 },
	{ "dwvwn",SF_FORMAT_DWVW_N },
	{ "dpcm8",SF_FORMAT_DPCM_8 },
	{ "dpcm16",SF_FORMAT_DPCM_16 },
	{ "vorbis",SF_FORMAT_VORBIS },
	{ "alac16",SF_FORMAT_ALAC_16 },
	{ "alac20",SF_FORMAT_ALAC_20 },
	{ "alac24",SF_FORMAT_ALAC_24 },
	{ "alac32",SF_FORMAT_ALAC_32 }
};

// map of default (ie sensible) subformats (usually 16-bit PCM)
const std::map<std::string, std::string> defaultSubFormats = {
	{"aiff","16"},
	{"au","16"},
	{"avr","16"},
	{"caf","16"},
	{"flac","16"},
	{"htk","16"},
	{"iff","16"},
	{"mat","16"},
	{"mpc","16"},
	{"oga","vorbis"},
	{"paf","16"},
	{"pvf","16"},
	{"raw","16"},
	{"rf64","16"},
	{"sd2","16"},
	{"sds","16"},
	{"sf","16"},
	{"voc","16"},
	{"w64","16"},
	{"wav","16"},
	{"wve","alaw"},
	{"xi","dpcm16"}
};

typedef enum {
	relaxed,
	normal,
	steep
} LPFMode;

// structure for holding all the parameters required for a conversion job:
struct conversionInfo
{
	std::string InputFilename;
	std::string OutputFilename;
	unsigned int OutputSampleRate;
	double Limit;
	bool bNormalize;
	int OutputFormat;
	bool bDither;
	double DitherAmount;
	int ditherProfileID;
	bool bAutoBlankingEnabled;
	bool bMinPhase;
	bool bSetFlacCompression;
	int flacCompressionLevel;
	bool bSetVorbisQuality;
	double vorbisQuality;
	bool disableClippingProtection;
	LPFMode lpfMode;
	bool bUseSeed;
	int seed;
	bool dsfInput;
};

bool determineBestBitFormat(std::string & BitFormat, const std::string & inFilename, const std::string & outFilename);
int determineOutputFormat(const std::string & outFileExt, const std::string & bitFormat);
void listSubFormats(const std::string & f);
int gcd(int a, int b);
Fraction GetSimplifiedFraction(int InputSampleRate, int OutputSampleRate);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, std::string & Parameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, unsigned int & nParameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, int & nParameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, double & Parameter);
bool findCmdlineOption(char ** begin, char ** end, const std::string & option);
template<typename FloatType> bool Convert(const conversionInfo& ci);
template<typename FloatType> bool dsfConvert(const conversionInfo & ci);

// Timer macros:
#define START_TIMER() LARGE_INTEGER starttime,finishtime,elapsed,frequency,timetaken; \
	QueryPerformanceFrequency(&frequency); \
	QueryPerformanceCounter(&starttime)

#define STOP_TIMER() QueryPerformanceCounter(&finishtime); \
	elapsed.QuadPart=finishtime.QuadPart-starttime.QuadPart; \
	timetaken.QuadPart=((1000*elapsed.QuadPart)/frequency.QuadPart); \
	std::cout << "Time=" << static_cast<long>(timetaken.QuadPart) << " ms" << std::endl

#endif // !RESAMPLER_H
