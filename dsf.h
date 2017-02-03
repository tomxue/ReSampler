#ifndef DSF_H_
#define DSF_H_

#include <cassert>
#include <cstdint>
#include <string>
#include <fstream>

#pragma pack(push, r1, 1)
typedef struct {
	uint32_t header;	// expected: "DSD "
	uint64_t length;	// expected: 28
	uint64_t filesize;
	uint64_t metadataPtr;
} DsfDSDChunk;

typedef enum {
	mono = 1,
	stereo,
	ch3,
	quad,
	ch4,
	ch5,
	ch51
} DsfChannelType;

typedef struct {
	uint32_t header;	// expected: "fmt "
	uint64_t length;	// expected: 52
	uint32_t version;	// expected: 1
	uint32_t formatID;	// expected: 0
	uint32_t channelType;
	uint32_t numChannels;
	uint32_t sampleRate;	// expected: 2822400 or 5644800
	uint32_t bitDepth;		// expected: 1 or 8
	uint64_t numSamples;
	uint32_t blockSize;	// expected: 4096
	uint32_t reserved;	// expected: zero
} DsfFmtChunk;

typedef struct {
	uint32_t header;	// expected: "data"
	uint64_t length;	// expected: 12 + sample data length
} DsfDataChunk;
#pragma pack(pop, r1)

typedef enum {
	dsf_read,
	dsf_write
} OpenMode;


// DsfFile interface:

class DsfFile 
{
public:
	// Construction / destruction
	DsfFile(const std::string& path, OpenMode mode = dsf_read) : path(path), mode(mode)
	{
		checkSizes();
		file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

		switch (mode) {
		case dsf_read:
			try {
				file.open(path, std::ios::in | std::ios::binary);
				readHeaders();
				for (int n = 0; n < 6; ++n) {
					channelBuffer[n] = new uint8_t[blockSize];
				}
				bufferIndex = blockSize; // initial state: empty (zero == full)
				currentBit = 0;
				currentChannel = 0;
			}
			catch (std::ios_base::failure& e) {
				std::cerr << e.what() << '\n';
			}

			break;
		case dsf_write:
			break;
		}
	};

	~DsfFile() {
		file.close();
		for (int n = 0; n < 6; ++n) {
			delete[] channelBuffer[n];
		}
	}

	// API:

	unsigned int channels() {
		return numChannels;
	};

	unsigned int sampleRate() {
		return _sampleRate;
	};

	int format() {
		return 0;
	};

	uint64_t frames() {
		return numFrames;
	};

	uint64_t samples() {
		return numSamples;
	};

	uint64_t read(float* buffer, uint64_t count) {

		// Interleaving in a dsf file is done at the block level.
		// This means that we read from file like this:
		// blockSize bytes of channel 0, blockSize bytes of channel 1 ... blockSize bytes of channel n 
		
		// However, caller expects interleaving to be done at the _sample_ level 

		uint64_t samplesRead = 0i64;

		for (uint64_t i = 0; i < count; ++i) {

			if (bufferIndex == blockSize) { // end of buffer ; fetch more data from file
				if (readBlocks() == 0) { // no more data
					break;
				}
				bufferIndex = 0;
			}

			buffer[i] = (channelBuffer[currentChannel][bufferIndex] & (1 << currentBit)) ? 1.0 : 0.0;
			++samplesRead;

			// cycle through channels, then bits, then bufferIndex

			if (++currentChannel == numChannels) {
				currentChannel = 0;
				if (++currentBit == 8) {
					currentBit = 0;
					++bufferIndex;
				}
			}
		}
		return samplesRead;
	};

	// testRead() : reads the entire file 
	// and confirms number of samples read equals number of samples expected:

	void testRead() {
		float sampleBuffer[8192];
		uint64_t totalSamplesRead = 0i64;
		uint64_t samplesRead = 0i64;

		while ((samplesRead = read(sampleBuffer, 8192)) != 0) {
			totalSamplesRead += samplesRead;
		}
		std::cout << "samples expected: " << numSamples << std::endl;
		std::cout << "total samples retrieved: " << totalSamplesRead << std::endl;
	}

private:
	DsfDSDChunk dsfDSDChunk;
	DsfFmtChunk dsfFmtChunk;
	DsfDataChunk dsfDataChunk;
	std::string path;
	OpenMode mode;
	std::fstream file;

	uint32_t blockSize;
	uint32_t numChannels;
	uint32_t _sampleRate;
	uint64_t numSamples;
	uint64_t numFrames;
	uint8_t* channelBuffer[6];
	uint64_t bufferIndex;
	uint32_t currentChannel;
	uint32_t currentBit;
	uint64_t endOfData;
	
	void checkSizes() {
		assert(sizeof(dsfDSDChunk) == 28);
		assert(sizeof(dsfFmtChunk) == 52);
		assert(sizeof(dsfDataChunk) == 12);
	}

	void readHeaders() {
		file.read((char*)&dsfDSDChunk, sizeof(dsfDSDChunk));
		file.read((char*)&dsfFmtChunk, sizeof(dsfFmtChunk));
		file.read((char*)&dsfDataChunk, sizeof(dsfDataChunk));

		blockSize = dsfFmtChunk.blockSize;
		numChannels = dsfFmtChunk.numChannels;
		_sampleRate = dsfFmtChunk.sampleRate;
		numFrames = dsfFmtChunk.numSamples;
		numSamples = numFrames * numChannels;
		endOfData = dsfDSDChunk.length + dsfFmtChunk.length + dsfDataChunk.length;

		assert( // metadata tag either non-existent or at end of data
			(dsfDSDChunk.metadataPtr == 0) ||
			(dsfDSDChunk.metadataPtr == endOfData)
			);
	}

	uint32_t readBlocks() {
		if (file.tellg() >= endOfData)
			return 0;

		for (int ch = 0; ch < numChannels; ++ch) {
			file.read((char*)channelBuffer[ch], blockSize);
		}
		return blockSize;
	}

};


/*

// how it's used in main():
unsigned int nChannels = infile.channels();
unsigned int InputSampleRate = infile.samplerate();
int InputFileFormat = infile.format();
sf_count_t InputSampleCount = infile.frames() * nChannels;

count = infile.read(inbuffer, BufferSize);

infile.seek(0i64, SEEK_SET); // rewind back to start of file

*/

#endif

