#ifndef IQDEMODULATOR_H
#define IQDEMODULATOR_H


#include <string>
#include <cstdint>
#include <memory>
#include <vector>
#include <iostream>

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

namespace  ReSampler {

enum ModulationType
{
	ModulationTypeNone,
	NFM,
	AM,
	LSB,
	USB,
	WFM,
	DSB,
	CW
};


class IQFile
{

public:

	IQFile(const std::string& fileName) : sndfileHandle(new SndfileHandle(fileName))
	{
	}

	IQFile(const std::string& fileName, int infileMode, int infileFormat, int infileChannels, int infileRate) : sndfileHandle(new SndfileHandle(fileName, infileMode, infileFormat, infileChannels, infileRate))
	{
	}

	bool error() {
		return sndfileHandle == nullptr || sndfileHandle->error();
	}

	int channels() {
		return (modulationType == ModulationType::WFM) ? 2 : 1;  // WFM is the only modulation type which can produce stereo
	}

	int samplerate() {
		return sndfileHandle == nullptr ? 0 : sndfileHandle->samplerate();
	}

	int64_t frames() {
		return sndfileHandle == nullptr ? 0LL : sndfileHandle->frames();
	}

	int format() {
		return sndfileHandle == nullptr ? 0 : sndfileHandle->format();
	}

	template<typename FloatType>
	int64_t read(FloatType* inbuffer, int64_t count) {

		if(sndfileHandle == nullptr) {
			std::cout << "nullptr" << std::endl;
			return 0LL;
		}

		if(wavBuffer.size() < count * 2) {
			wavBuffer.resize(count * 2);
		}

		int64_t samplesRead = sndfileHandle->read(wavBuffer.data(), count * 2);

		int64_t j = 0;
		for(int64_t i = 0; i < samplesRead; i += 2) {
			inbuffer[j++] = demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1));
		}

		return j;
	}

	sf_count_t seek(int64_t frames, int whence) {
		if(sndfileHandle == nullptr) {
			return 0LL;
		}

		return sndfileHandle->seek(frames, whence);
	}

// getters
	ModulationType getModulationType() const
	{
		return modulationType;
	}

// setters
	void setModulationType(const ModulationType &value)
	{
		modulationType = value;
	}

private:
	template<typename FloatType>
	FloatType demodulateFM(FloatType i, FloatType q)
	{
		// this is actually quite simple, thanks to some clever calculus tricks.
		// see https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/

		i2 = i1;
		i1 = i0;
        i0 = i;
        q2 = q1;
        q1 = q0;
        q0 = q;

		double gain = 1.0 / (0.001 + i * i + q * q); // todo: calcuate correct bias amount to avoid div / 0
		return gain * (((q0 - q2) * i1) - ((i0 - i2) * q1));
    }

private:
	// resources
	std::unique_ptr<SndfileHandle> sndfileHandle;
	std::vector<double> wavBuffer;

	// properties
	ModulationType modulationType{ModulationType::ModulationTypeNone};

	// state
	double i0{0.0};
	double i1{0.0};
	double i2{0.0};
	double q0{0.0};
	double q1{0.0};
	double q2{0.0};

};


} // namespace  ReSampler

#endif // IQDEMODULATOR_H
