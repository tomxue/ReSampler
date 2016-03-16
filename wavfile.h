#ifndef WAVFILE_H
#define WAVFILE_H 1

///////////////////////////////////////////////
// Wave File I/O library (C) Judd Niemann 2016
//
//////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#define WAVE_FORMAT_PCM 0x0001
#define WAVE_FORMAT_IEEE_FLOAT 0x0003
#define WAVE_FORMAT_ALAW 0x0006					// 8 - bit ITU - T G.711 A - law
#define WAVE_FORMAT_MULAW 0x0007				// 8 - bit ITU - T G.711 µ - law
#define WAVE_FORMAT_EXTENSIBLE 0xFFFE			// Determined by SubFormat

typedef struct WavHeader
{
	unsigned __int32 ChunkID;		// "RIFF" (BIG ENDIAN)
	unsigned __int32 ChunkSize;		// 4+ size of all subsequent sub Chunks (iow, size of remainder of file (filesize-8))
	unsigned __int32 Format;		// "WAVE" (BIG ENDIAN)

	unsigned __int32 fmtChunkID;	// "fmt " (BIG ENDIAN)
	unsigned __int32 fmtChunkSize; // 16,18, or 40
	
	unsigned __int16 AudioFormat;
	unsigned __int16 NumChannels;
	unsigned __int32 SampleRate;
	unsigned __int32 ByteRate;
	unsigned __int16 BlockAlign;
	unsigned __int16 BitsPerSample;
	
	// only used if (fmtChunkSize == 18 || fmtChunkSize == 40) :
	unsigned __int16 cbSize; 
	
	// only used if (fmtChunkSize == 40 && cbSize==22) :
	unsigned __int16 ValidBitsPerSample;
	unsigned __int32 ChannelMask;
	unsigned char SubFormat[16];

	// Only used if (AudioFormat != WAVE_FORMAT_PCM) :
	unsigned __int32 factChunkID;	// "fact" (BIG ENDIAN)
	unsigned __int32 factChunkSize;
	unsigned __int32 SampleLength;
	
	unsigned __int32 dataChunkID;	// "data" (BIG ENDIAN)
	unsigned __int32 dataChunkSize;

} WavHeader;

template<typename SampleType> class InputWavFile {
public:
	
	bool Open(const std::string& filename) {
		// To-do: Open file
		infile.open(filename, std::ios::binary | std::ios::in);
		memset(static_cast<void*>(&Header), 0, sizeof(Header));
		return readHeader();
	}
	unsigned __int64 Read(SampleType* buffer, unsigned int nSamples) {}
	void Close() {
		infile.close();
	}
private:
	std::ifstream infile;
	WavHeader Header;
	bool readHeader(){
		auto savedPosition = infile.tellg();
		
		// expected: master chunk:
		infile.read(reinterpret_cast<char*>(&Header.ChunkID), 4);
		infile.read(reinterpret_cast<char*>(&Header.ChunkSize), 4);
		infile.read(reinterpret_cast<char*>(&Header.Format), 4);
		
		if (Header.ChunkID != 0x46464952	/* "FFIR" */ 
			|| Header.Format != 0x45564157	/* "EVAW" */ ) { // Not .wav file
			infile.close();
			return false;
		}

		// expected: "fmt " sub-chunk:
		infile.read(reinterpret_cast<char*>(&Header.fmtChunkID),4);
		if (Header.fmtChunkID != 0x20746d66 /* " tmf" */) {
			infile.close();
			return false;
		}
		else {
			infile.read(reinterpret_cast<char*>(&Header.fmtChunkSize), 4); // 16,18, or 40
			infile.read(reinterpret_cast<char*>(&Header.AudioFormat), 2);
			infile.read(reinterpret_cast<char*>(&Header.NumChannels), 2);
			infile.read(reinterpret_cast<char*>(&Header.SampleRate), 4);
			infile.read(reinterpret_cast<char*>(&Header.ByteRate), 4);
			infile.read(reinterpret_cast<char*>(&Header.BlockAlign), 2);
			infile.read(reinterpret_cast<char*>(&Header.BitsPerSample), 2);

			if (Header.fmtChunkSize == 18 || Header.fmtChunkSize == 40){
				unsigned __int16 cbSize;
				if (Header.fmtChunkSize == 40 && Header.cbSize == 22) {
					unsigned __int16 ValidBitsPerSample;
					unsigned __int32 ChannelMask;
					unsigned char SubFormat[16];
				}
			}
		}
		
		switch (Header.AudioFormat) {
		case WAVE_FORMAT_PCM:
			break;
		case WAVE_FORMAT_IEEE_FLOAT:
		case WAVE_FORMAT_ALAW:
		case WAVE_FORMAT_MULAW:
		case WAVE_FORMAT_EXTENSIBLE:
		default:
			// for non-PCM formats, expected: "fact" SubChunk
			savedPosition = infile.tellg();
			infile.read(reinterpret_cast<char*>(&Header.factChunkID), 4);
			if (Header.factChunkID == 0x74636166 /* "tcaf" */) { 
				infile.read(reinterpret_cast<char*>(&Header.factChunkSize), 4);
				infile.read(reinterpret_cast<char*>(&Header.SampleLength), 4);
			}
			else { // not a "fact" SubChunk (non-fatal)
				infile.seekg(savedPosition); // back up to start of SubChunk ID
			}
		}

		// skim through remaining (unsupported) sub chunks:
		do { 
			__int32 ChunkID;
			__int32 ChunkSize;
			savedPosition = infile.tellg();
			
			infile.read(reinterpret_cast<char*>(&ChunkID), 4);
			if (ChunkID == 0x61746164 /* "atad" */) { // Found "data" sub chunk
				infile.seekg(savedPosition);
				break;
			}
			infile.read(reinterpret_cast<char*>(&ChunkSize), 4);
			infile.seekg(ChunkSize - 4 + infile.tellg());
		} while (!infile.eof());

		// expected: "data" subChunk
		infile.read(reinterpret_cast<char*>(&Header.dataChunkID), 4);
		if (Header.dataChunkID != 0x61746164 /* "atad" */ ) {
			infile.close();
			return false;
		}
		infile.read(reinterpret_cast<char*>(&Header.dataChunkSize), 4);
		return true;
	}
};


template<typename SampleType> class OutputWavFile {
public:
	bool Open(const std::string& filename) {}
	unsigned __int64 Write(const SampleType* buffer, unsigned int nSamples) {}
	bool Close() {}
};

#endif

