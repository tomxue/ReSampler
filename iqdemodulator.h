#ifndef IQDEMODULATOR_H
#define IQDEMODULATOR_H

#include <string>
#include <cstdint>
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>
#include <map>

#include <sndfile.h>
#include <sndfile.hh>

namespace  ReSampler {

enum ModulationType
{
	ModulationTypeNone = 0,
	NFM,
	AM,
	LSB,
	USB,
	WFM,
	DSB,
	CW
};

static const std::map<std::string, ModulationType> modulationTypeMap
{
	{"NONE", ModulationType::ModulationTypeNone},
	{"NFM", ModulationType::NFM},
	{"AM", ModulationType::AM},
	{"LSB", ModulationType::LSB},
	{"USB", ModulationType::USB},
	{"WFM", ModulationType::WFM},
	{"DSB", ModulationType::DSB},
	{"CW", ModulationType::CW}
};

class IQFile
{
public:

	IQFile(const std::string& fileName) : sndfileHandle(new SndfileHandle(fileName))
	{
	}

	IQFile(const std::string& fileName, int infileMode, int infileFormat, int infileChannels, int infileRate) :
		sndfileHandle(new SndfileHandle(fileName, infileMode, infileFormat & 0xFFFF00FF, infileChannels, infileRate))
	{
        // Extract modulation type from 2nd-last byte of file format code.
        // (Note: libsndfile has this for the subformat mask:
        // SF_FORMAT_SUBMASK = 0x0000FFFF
        // So far, only the least-significant byte has been used. ie: 0x000000FF.
        // If they ever add more formats in the future which use the upper byte,
        // then this strategy may need reevaluation ...)

		modulationType = static_cast<ModulationType>((infileFormat & 0x0000FF00) >>  8);
	}

	bool error() {
		return sndfileHandle == nullptr || sndfileHandle->error();
	}

	int channels() {
		return 1;  // I & Q inputs always get demodulated into a single channel
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
			switch(modulationType) {
			case AM:
				inbuffer[j++] = demodulateAM(wavBuffer.at(i), wavBuffer.at(i + 1));
				break;
            case LSB:
            case USB:
                // SSB : just copy I-component
                inbuffer[j++] = wavBuffer.at(i);
                break;
			default:
                // Narrowband FM
				inbuffer[j++] = demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1));
			}
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
		static constexpr double maxGain = 60.0;
		static const double c = std::pow(10.0, -(maxGain / 20.0));

		// this is actually quite simple, thanks to some clever calculus tricks.
		// see https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/

		i2 = i1;
		i1 = i0;
        i0 = i;
        q2 = q1;
        q1 = q0;
        q0 = q;

		double gain = 1.0 / (c + i * i + q * q);
		return gain * (((q0 - q2) * i1) - ((i0 - i2) * q1));
    }

	template<typename FloatType>
	FloatType demodulateAM(FloatType i, FloatType q)
	{
		static constexpr FloatType scale = 0.7071; // << 1/sqrt(2)
		return scale * std::sqrt(i * i + q * q);
	}

private:
	// resources
	std::unique_ptr<SndfileHandle> sndfileHandle;
	std::vector<double> wavBuffer;

	// properties
	ModulationType modulationType{ModulationType::NFM};

	// registers used for demodulating FM
	double i0{0.0};
	double i1{0.0};
	double i2{0.0};
	double q0{0.0};
	double q1{0.0};
	double q2{0.0};

};

} // namespace  ReSampler

#endif // IQDEMODULATOR_H
