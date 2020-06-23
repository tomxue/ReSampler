#ifndef IQDEMODULATOR_H
#define IQDEMODULATOR_H


#include <string>
#include <cstdint>
#include <sndfile.h>
#include <sndfile.hh>

/* Note: type 'FileReader' MUST implement the following methods:
constructor(const std::string& fileName)
constructor(const std::string& fileName, int infileMode, int infileFormat, int infileChannels, int infileRate)
bool error() // or int error()
unsigned int channels()
unsigned int samplerate()
uint64_t frames()
int format()
read(inbuffer, count)
seek(position, whence)
*/


template<typename FloatType>
class IQFile
{

public:

	IQFile(const std::string& fileName) {
		(void)fileName;
	}

	IQFile(const std::string& fileName, int infileMode, int infileFormat, int infileChannels, int infileRate) {
		(void)fileName;
		(void)infileMode;
		(void)infileFormat;
		(void)infileChannels;
		(void)infileRate;
	}

	bool error() {
		return false;
	}

	int channels() {
		return 0;
	}

	int samplerate() {
		return 0;
	}

	int64_t frames() {
		return 0LL;
	}

	int format() {
		return 0;
	}

	int64_t read(FloatType* inbuffer, int64_t count) {
		(void)inbuffer;
		(void)count;
		return 0LL;
	}

	sf_count_t seek(int64_t frames, int whence) {
		(void)frames;
		(void)whence;
		return 0LL;
	}

private:
    FloatType demodulateFM(FloatType i, FloatType q)
    {
        i2 = i1;
        i1 = i0;
        i0 = i;
        q2 = q1;
        q1 = q0;
        q0 = q;

        return ((q0 - q2) * i1) - ((i0 - i2) * q1);
    }

private:
    FloatType i0{0.0};
    FloatType i1{0.0};
    FloatType i2{0.0};
    FloatType q0{0.0};
    FloatType q1{0.0};
    FloatType q2{0.0};

};

#endif // IQDEMODULATOR_H
