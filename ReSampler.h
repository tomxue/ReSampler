/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef RESAMPLER_H
#define RESAMPLER_H 1

#include "sndfile.h"
#include "sndfile.hh"

#include <string>
#include <map>

struct ConversionInfo;

const std::string strVersion("1.3.7");
const std::string strUsage("usage: ReSampler -i <inputfile> [-o <outputfile>] -r <samplerate> [-b <bitformat>] [-n [<normalization factor>]]\n");
const std::string strExtraOptions(
	"--help\n"
	"--version\n"
	"--compiler\n"
	"--sndfile-version\n"
	"--listsubformats <ext>\n"
	"--showDitherProfiles\n"
	"--gain [<amount>]\n"
	"--doubleprecision\n"
	"--dither [<amount>] [--autoblank] [--ns [<ID>]] [--flat-tpdf] [--seed [<num>]]\n"
	"--noDelayTrim\n"
	"--minphase\n"
	"--flacCompression <compressionlevel>\n"
	"--vorbisQuality <quality>\n"
	"--noClippingProtection\n"
	"--relaxedLPF\n"
	"--steepLPF\n"
	"--lpf-cutoff <percentage> [--lpf-transition <percentage>]\n"
	"--mt\n"
	"--rf64\n"
	"--noPeakChunk\n"
	"--noMetadata\n"
	"--maxStages\n"
	"--showStages\n"
);
const double clippingTrim = 1.0 - (1.0 / (1 << 24));

#define BUFFERSIZE 32768 // buffer size for file reads
#define MAXCHANNELS 64

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

#define MAX_CART_TAG_TEXT_SIZE 16384
typedef SF_CART_INFO_VAR(MAX_CART_TAG_TEXT_SIZE) LargeSFCartInfo;

typedef struct
{
	std::string title;
	std::string copyright;
	std::string software;
	std::string artist;
	std::string comment;
	std::string date;
	std::string album;
	std::string license;
	std::string trackNumber;
	std::string genre;

	// The following is only relevant for bext chunks in Broadcast Wave files:
	bool has_bext_fields;
	SF_BROADCAST_INFO broadcastInfo;

	// The following is only relevant for cart chunks:
	bool has_cart_chunk;
	LargeSFCartInfo cartInfo;

} MetaData;

bool checkSSE2();
bool checkAVX();
bool showBuildVersion();
bool parseParameters(ConversionInfo & ci, bool & bBadParams, int argc, char * argv[]);
bool determineBestBitFormat(std::string & BitFormat, const std::string & inFilename, const std::string & outFilename);
int determineOutputFormat(const std::string & outFileExt, const std::string & bitFormat);
void listSubFormats(const std::string & f);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, std::string & Parameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, unsigned int & nParameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, int & nParameter);
void getCmdlineParam(char ** begin, char ** end, const std::string & OptionName, double & Parameter);
bool findCmdlineOption(char ** begin, char ** end, const std::string & option);
template<typename FileReader, typename FloatType> bool convert(ConversionInfo & ci, bool peakDetection = true);
int getDefaultNoiseShape(int sampleRate);
void showDitherProfiles();
int getSfBytesPerSample(int format);
bool checkWarnOutputSize(uint64_t inputSamples, int bytesPerSample, int numerator, int denominator);
std::string fmtNumberWithCommas(uint64_t n);
bool getMetaData(MetaData& metadata, SndfileHandle& infile);
bool setMetaData(const MetaData& metadata, SndfileHandle& outfile);
void showCompiler();

#endif // !RESAMPLER_H