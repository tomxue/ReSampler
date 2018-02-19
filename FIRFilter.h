/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef FIRFFILTER_H_
#define FIRFFILTER_H_

// FIRFilter.h : simple FIR filter implementation by J.Niemann

// #define USE_SIMD_FOR_DOUBLES

#include <typeinfo>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <cassert>
#include <xmmintrin.h>

#include <fftw3.h>
#include "alignedmalloc.h"

#define FILTERSIZE_LIMIT 131071
#define FILTERSIZE_BASE 103

#define SSE_ALIGNMENT_SIZE 16

#ifdef USE_AVX
#include "FIRFilterAVX.h"

#else
#if (defined(_M_X64) || defined(__x86_64__) || defined(USE_SSE2)) // All x64 CPUs have SSE2 instructions, but some older 32-bit CPUs do not. 
	#define USE_SIMD 1 // Vectorise main loop in FIRFilter::get() by using SSE2 SIMD instrinsics 
#endif

#if defined (__MINGW64__) || defined (__MINGW32__) || defined (__GNUC__)
#ifdef USE_QUADMATH
#include <quadmath.h>
#define FIR_QUAD_PRECISION
#endif
#endif

template <typename FloatType>
class FIRFilter {

public:

	// constructor:
	FIRFilter(const FloatType* taps, size_t size) :
		size(size), Signal(nullptr), CurrentIndex(size-1), LastPut(0),
		Kernel0(nullptr), Kernel1(nullptr), Kernel2(nullptr), Kernel3(nullptr)
	{

		// allocate buffers:
		allocateBuffers();
		assertAlignment();

		// initialize filter kernel and signal buffers
		for (unsigned int i = 0; i < size; ++i) {
			Kernel0[i] = taps[i];
			Signal[i] = 0.0;
			Signal[i + size] = 0.0;
		}

		// Populate additional kernel Phases:
		memcpy(1 + Kernel1, Kernel0, (size - 1) * sizeof(FloatType));
		Kernel1[0] = Kernel0[size - 1];
		memcpy(1 + Kernel2, Kernel1, (size - 1) * sizeof(FloatType));
		Kernel2[0] = Kernel1[size - 1];
		memcpy(1 + Kernel3, Kernel2, (size - 1) * sizeof(FloatType));
		Kernel3[0] = Kernel2[size - 1];
	}

	// deconstructor:
	~FIRFilter() {
		freeBuffers();
	}

	// copy constructor:
	FIRFilter(const FIRFilter& other) : size(other.size), CurrentIndex(other.CurrentIndex), LastPut(other.LastPut)
	{
		allocateBuffers();
		assertAlignment();
		copyBuffers(other);
	}

	// move constructor:
	FIRFilter(FIRFilter&& other) noexcept :
		size(other.size), Signal(other.Signal), CurrentIndex(other.CurrentIndex), LastPut(other.LastPut),
		Kernel0(other.Kernel0), Kernel1(other.Kernel1), Kernel2(other.Kernel2), Kernel3(other.Kernel3)
	{
		assertAlignment();
		other.Signal = nullptr;
		other.Kernel0 = nullptr;
		other.Kernel1 = nullptr;
		other.Kernel2 = nullptr;
		other.Kernel3 = nullptr;
	}
	
	// copy assignment:
	FIRFilter& operator= (const FIRFilter& other)
	{
		size = other.size;
		CurrentIndex = other.CurrentIndex;
		LastPut = other.LastPut;
		freeBuffers();
		allocateBuffers();
		assertAlignment();
		copyBuffers(other);
		return *this;
	}

	// move assignment:
	FIRFilter& operator= (FIRFilter&& other) noexcept
	{
		if(this!=&other) // prevent self-assignment
		{
			size = other.size;
			CurrentIndex = other.CurrentIndex;
			LastPut = other.LastPut;

			freeBuffers();
			
			Signal = other.Signal;
			Kernel0 = other.Kernel0;
			Kernel1 = other.Kernel1;
			Kernel2 = other.Kernel2;
			Kernel3 = other.Kernel3;
			assertAlignment();

			other.Signal = nullptr;
			other.Kernel0 = nullptr;
			other.Kernel1 = nullptr;
			other.Kernel2 = nullptr;
			other.Kernel3 = nullptr;
		}
		return *this;
	}

	bool operator== (const FIRFilter& other) const 
	{
		if (size != other.size)
			return false;
		
		for (int i = 0; i < size; i++) {
			if (Kernel0[i] != other.Kernel0[i])
				return false;
		}
		
		return true;
	}

	void reset() {
		// reset indexes:
		CurrentIndex = size - 1;
		LastPut = 0;

		// clear signal buffer
		for (unsigned int i = 0; i < size; ++i) {
			Signal[i] = 0.0;
			Signal[i + size] = 0.0;
		}

	}

	void put(FloatType value) { // Put signal in reverse order.
		Signal[CurrentIndex] = value;
		LastPut = CurrentIndex;
		if (CurrentIndex == 0) {
			CurrentIndex = size - 1; // Wrap
			memcpy(Signal + size, Signal, size*sizeof(FloatType)); // copy history to upper half of buffer
		}
		else
			--CurrentIndex;
	}

	void putZero() {
		Signal[CurrentIndex] = 0.0;
		if (CurrentIndex == 0) {
			CurrentIndex = size - 1; // Wrap
			memcpy(Signal + size, Signal, size*sizeof(FloatType)); // copy history to upper half of buffer
		}
		else
			--CurrentIndex;
	}

	FloatType get() {

#ifndef USE_SIMD
		FloatType output = 0.0;
		int index = CurrentIndex;
		for (int i = 0; i < size; ++i) {
			output += Signal[index] * Kernel0[i];
			index++;
		}
		return output;
#else
		// SIMD implementation: This only works with floats (doubles need specialisation)

		FloatType output = 0.0;
		FloatType* Kernel = Kernel0;
		int Index = (CurrentIndex >> 2) << 2; // make multiple-of-four
		int Phase = CurrentIndex & 3;
		
		// Part1 : Head
		// select proper Kernel phase and calculate first Block of 4:
		switch (Phase) {
		case 0:
			Kernel = Kernel0;
			// signal already aligned and ready to use
			output = Kernel[0] * Signal[Index] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
			break;
		case 1:
			Kernel = Kernel1;
			// signal starts at +1 : load first value from history (ie upper half of buffer)
			output = Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
			break;
		case 2:
			Kernel = Kernel2;
			// signal starts at +2 : load first and second values from history (ie upper half of buffer)
			output = Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
			break;
		case 3: 
			Kernel = Kernel3;
			// signal starts at +3 : load 1st, 2nd, 3rd values from history (ie upper half of buffer)
			output = Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + 3];
			break;
		}
		Index += 4;

		// Part 2: Body
		alignas(SSE_ALIGNMENT_SIZE) __m128 signal;	// SIMD Vector Registers for calculation
		alignas(SSE_ALIGNMENT_SIZE) __m128 kernel;
		alignas(SSE_ALIGNMENT_SIZE) __m128 product;
		alignas(SSE_ALIGNMENT_SIZE) __m128 accumulator = _mm_setzero_ps();

		for (int i = 4; i < (size >> 2) << 2; i += 4) {
			signal = _mm_load_ps(Signal + Index);
			kernel = _mm_load_ps(Kernel + i);
			product = _mm_mul_ps(signal, kernel);
			accumulator = _mm_add_ps(product, accumulator);
			Index += 4;
		}
		
		// http://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
		__m128 a   = _mm_shuffle_ps(
			accumulator, 
			accumulator,                 // accumulator = [D     C     | B     A    ]
			_MM_SHUFFLE(2, 3, 0, 1));                  // [C     D     | A     B    ]
		__m128 b   = _mm_add_ps(accumulator, a);       // [D+C   C+D   | B+A   A+B  ]
		a          = _mm_movehl_ps(a, b);              // [C     D     | D+C   C+D  ]
		b          = _mm_add_ss(a, b);                 // [C     D     | D+C A+B+C+D]
		output    += _mm_cvtss_f32(b);                 // A+B+C+D
	
		// Part 3: Tail
		for (int j = (size >> 2) << 2; j < size; ++j) {
			output += Signal[Index] * Kernel[j];
			++Index;
		}

		return output;

#endif // !USE_SIMD
	}

	FloatType lazyGet(int L) {	// Skips stuffed-zeros introduced by interpolation, by only calculating every Lth sample from LastPut
		FloatType output = 0.0;
		int Offset = LastPut - CurrentIndex;
		if (Offset < 0) { // Wrap condition
			Offset += size;
		}
	
		for (int i = Offset; i < size; i+=L) {
			output += Signal[i+ CurrentIndex] * Kernel0[i];
		}
		return output;
	}

private:
	size_t size;
	FloatType* Signal; // Double-length signal buffer, to facilitate fast emulation of a circular buffer
	int CurrentIndex;
	int LastPut;

	// Polyphase Filter Kernel table:
	FloatType* Kernel0;
	FloatType* Kernel1;
	FloatType* Kernel2;
	FloatType* Kernel3;

	void allocateBuffers()
	{
		Signal = static_cast<FloatType*>(aligned_malloc(2 * size * sizeof(FloatType), SSE_ALIGNMENT_SIZE));
		Kernel0 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), SSE_ALIGNMENT_SIZE));
		Kernel1 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), SSE_ALIGNMENT_SIZE));
		Kernel2 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), SSE_ALIGNMENT_SIZE));
		Kernel3 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), SSE_ALIGNMENT_SIZE));
	}

	void copyBuffers(const FIRFilter& other)
	{
		memcpy(Signal, other.Signal, 2 * size * sizeof(FloatType));
		memcpy(Kernel0, other.Kernel0, size * sizeof(FloatType));
		memcpy(Kernel1, other.Kernel1, size * sizeof(FloatType));
		memcpy(Kernel2, other.Kernel2, size * sizeof(FloatType));
		memcpy(Kernel3, other.Kernel3, size * sizeof(FloatType));
	}

	void freeBuffers()
	{
		aligned_free(Signal);
		aligned_free(Kernel0);
		aligned_free(Kernel1);
		aligned_free(Kernel2);
		aligned_free(Kernel3);
	}
	
	// assertAlignment() : asserts that all private data buffers are aligned on expected boundaries
	void assertAlignment()
	{
		const std::uintptr_t alignment = SSE_ALIGNMENT_SIZE;
		assert(reinterpret_cast<std::uintptr_t>(Signal) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(Kernel0) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(Kernel1) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(Kernel2) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(Kernel3) % alignment == 0);
	}
};

#ifdef USE_SIMD

#ifndef USE_SIMD_FOR_DOUBLES

// Specialization for doubles:
template <>
double FIRFilter<double>::get() {
	double output = 0.0;
	int index = CurrentIndex;
	for (int i = 0; i < size; ++i) {
		output += Signal[index] * Kernel0[i];
		index++;
	}
	return output;
}

#else 
 
// actual SIMD implementations for doubles (not worth the effort - no faster than than naive):

template <>
double FIRFilter<double>::get() {

	// SIMD implementation
	// Processes two doubles at a time.

	double output = 0.0;
	double* Kernel;
	int Index = (CurrentIndex >> 1) << 1; // make multiple-of-two
	int Phase = CurrentIndex & 1;

	// Part1 : Head
	// select proper Kernel phase and calculate first Block of 2:
	switch (Phase) {
	case 0:
		Kernel = Kernel0;
		// signal already aligned and ready to use
		output = Kernel[0] * Signal[Index] + Kernel[1] * Signal[Index + 1];
		break;
	case 1:
		Kernel = Kernel1;
		// signal starts at +1 : load first value from history (ie upper half of buffer)
		output = Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + 1];
		break;
	}
	Index += 2;

	// Part 2: Body
	alignas(SSE_ALIGNMENT_SIZE) __m128d signal;	// SIMD Vector Registers for calculation
	alignas(SSE_ALIGNMENT_SIZE) __m128d kernel;
	alignas(SSE_ALIGNMENT_SIZE) __m128d product;
	alignas(SSE_ALIGNMENT_SIZE) __m128d accumulator = _mm_setzero_pd();

	for (int i = 2; i < (size >> 1) << 1; i += 2) {
		signal = _mm_load_pd(Signal + Index);
		kernel = _mm_load_pd(Kernel + i);
		product = _mm_mul_pd(signal, kernel);
		accumulator = _mm_add_pd(product, accumulator);
		Index += 2;
	}

	output += accumulator.m128d_f64[0] + accumulator.m128d_f64[1];

	// Part 3: Tail
	for (int j = (size >> 1) << 1; j < size; ++j) {
		output += Signal[Index] * Kernel[j];
		++Index;
	}

	return output;
}
#endif // USE_SIMD_FOR_DOUBLES
#endif // USE_SIMD
#endif // !USE_AVX

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Functions beyond this point are for manipulating filter taps, and not for actually performing filtering -- //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// makeLPF() : generate low pass filter coefficients, using sinc function
template<typename FloatType> bool makeLPF(FloatType* filter, int Length, FloatType transitionFreq, FloatType sampleRate)
{
#ifdef FIR_QUAD_PRECISION

    std::cout << "calculating quad-precision filter coefficients ...\n";

    // use quads internally, regardless of FloatType
    __float128 ft = transitionFreq / sampleRate; // normalised transition frequency
    assert(ft < 0.5Q);
    int halfLength = Length / 2;
    __float128 halfM = 0.5Q * (Length - 1);

    if (Length & 1)
        filter[halfLength] = 2.0Q * ft; // if length is odd, avoid divide-by-zero at centre-tap

    for (int n = 0; n<halfLength; ++n) {
        __float128 sinc = sinq(2.0Q * M_PIq * ft * (n - halfM)) / (M_PIq * (n - halfM));	// sinc function
        filter[Length - n - 1] = filter[n] = sinc;	// exploit symmetry
    }

#else
	// use doubles internally, regardless of FloatType
	double ft = transitionFreq / sampleRate; // normalised transition frequency
	assert(ft < 0.5);
	int halfLength = Length / 2;
	double halfM = 0.5 * (Length - 1);

	if (Length & 1)
		filter[halfLength] = 2.0 * ft; // if length is odd, avoid divide-by-zero at centre-tap

	for (int n = 0; n<halfLength; ++n) {
		double sinc = sin(2.0 * M_PI * ft * (n - halfM)) / (M_PI * (n - halfM));	// sinc function
		filter[Length - n - 1] = filter[n] = sinc;	// exploit symmetry
	}
#endif

	return true;
}

// This function converts a requested sidelobe height (in dB) to a value for the Beta parameter used in a Kaiser window:
template<typename FloatType> FloatType calcKaiserBeta(FloatType dB) 
{
	if(dB<21.0)
	{
		return 0;
	}
	else if ((dB >= 21.0) && (dB <= 50.0)) {
		return 0.5842 * pow((dB - 21), 0.4) + 0.07886 * (dB - 21);
	}
	else if (dB>50.0) {
		return 0.1102 * (dB - 8.7);
	}
	else
	{
		return 0;
	}
}

// I0() : 0th-order Modified Bessel function of the first kind:
double I0(double z)
{
	double result = 0.0;
	double kfact = 1.0;
	for (int k = 0; k < 30; ++k) {
		if (k) kfact *= static_cast<double>(k);
		result += pow((pow(z, 2.0) / 4.0), k) / pow(kfact, 2.0);
	}
	// std::cout << "input: " << z << " output: " << result << std::endl;
	return result;
}

#ifdef FIR_QUAD_PRECISION
__float128 I0q(__float128 x)
{
    __float128 result = 0.0Q;
    __float128 kfact = 1.0Q;
    for (int k = 0; k < 60; ++k) {
        if (k)
            kfact *= k;
        result += powq((powq(x, 2.0Q) / 4.0Q), k) / powq(kfact, 2.0Q);
    }
    return result;
}
#endif

// applyKaiserWindow() - This function applies a Kaiser Window to an array of filter coefficients ("textbook" version):
template<typename FloatType> bool applyKaiserWindow(FloatType* filter, int Length, double Beta)
{
	// Note: sometimes, the Kaiser Window formula is defined in terms of Alpha (instead of Beta), 
	// in which case, Alpha def= Beta / pi

    if (Length < 1)
        return false;

#ifdef FIR_QUAD_PRECISION

    std::cout << "Applying quad-precision Kaiser Window ...\n";

    for (int n = 0; n < Length; ++n) {
        filter[n] *= I0q(Beta * sqrtq(1.0Q - powq((2.0Q * n / (Length - 1) - 1), 2.0Q)))
                     / I0q(Beta);
    }
#else
	for (int n = 0; n < Length; ++n) {
		filter[n] *= I0(Beta * sqrt(1.0 - pow((2.0 * n / (Length - 1) - 1), 2.0)))
			/ I0(Beta);
	}
#endif

	return true;
}

// applyKaiserWindow2() - applies a Kaiser Window to an array of filter coefficients (alternative formula):
template<typename FloatType> bool applyKaiserWindow2(FloatType* filter, int Length, double Beta)
 {
	 double A;	// use double internally, regardless of FloatType (speed not an issue here; no reason not to)
	 double maxA = 0; // for diagnostics
	 for (int n = 0; n < Length; ++n) {

		 // simplified Kaiser Window Equation:
		 A = (2.0 * Beta / Length) * sqrt(n*(Length - n - 1));
		 maxA = std::max(maxA, A);
		 filter[n] *= I0(A) / I0(Beta);
	 }
	
	//// diagnostic to check accuracy of I0():
	//	std::cout << "I0( " << maxA << " ) ==" << I0(maxA) << std::endl;
	
	return true;
}


// the following is a set of Complex-In, Complex-Out transforms used for constructing a minimum-Phase FIR:

// logV() : logarithm of a vector of Complex doubles
std::vector<std::complex<double>>
logV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return std::log(x); });
	return output;
}

// limitDynRangeV() : set a limit (-dB) on how quiet signal is allowed to be below the peak. 
// Guaranteed to never return zero.
std::vector<std::complex<double>>
limitDynRangeV(const std::vector<std::complex<double>>& input, double dynRangeDB) {
	double dynRangeLinear = pow(10, std::abs(dynRangeDB) / 20.0); // will give same result for positive or negative dB values.
	
	// find peak:
	double peak=0.0;
	for (auto &c : input) {
		peak = std::max(peak, std::abs(c));
	}
	
	// determine low threshold
	double lowThresh = peak / dynRangeLinear;	// a level which is dynRangeDB dB below peak
	std::complex<double> lastX = lowThresh;		// variable for storing last output value
	
	std::vector<std::complex<double>> output(input.size(), 0);

	std::transform(input.begin(), input.end(), output.begin(),
		[lowThresh, &lastX](std::complex<double> x) -> std::complex<double> {
		
		double level = std::abs(x);
		if (level < lowThresh) {
			if (level == 0.0) {		// when input is zero, we must somehow make the modulus of the output equal to lowThresh
				x = lastX;			// sticky output; use last output value instead of zero
			}
			else {
				x = (x / level) * lowThresh; // scale x such that |x| == lowThresh
				lastX = x;
			}
		}
		return x; } // ends lambda
	); // ends call to std::transform()

	return output;
}

// realV() : real parts of a vector of Complex doubles
std::vector<std::complex<double>>
realV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return x.real(); });
	return output;
}

// imagV() : imaginary parts of a vector of Complex doubles (answer placed in imaginary part of output):
std::vector<std::complex<double>>
imagV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return{ 0,x.imag() }; });
	return output;
}

// expV() : exp of a vector of Complex doubles
std::vector<std::complex<double>>
expV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return exp(x); });
	return output;
}

// fftV() : FFT of vector of Complex doubles
std::vector<std::complex<double>>
fftV(std::vector<std::complex<double>> input) {
	
	std::vector<std::complex<double>> output(input.size(), 0); // output vector
		
	// create, execute, destroy plan:
	fftw_plan p = fftw_plan_dft_1d(static_cast<int>(input.size()), 
		reinterpret_cast<fftw_complex*>(&input[0]), 
		reinterpret_cast<fftw_complex*>(&output[0]), 
		FFTW_FORWARD, 
		FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);
	
	return output;
}

// ifftV() : Inverse FFT of vector of Complex doubles
std::vector<std::complex<double>>
ifftV(std::vector<std::complex<double>> input) {

	std::vector<std::complex<double>> output(input.size(), 0); // output vector

	// create, execute, destroy plan:
	fftw_plan p = fftw_plan_dft_1d(static_cast<int>(input.size()),
		reinterpret_cast<fftw_complex*>(&input[0]),
		reinterpret_cast<fftw_complex*>(&output[0]),
		FFTW_BACKWARD,
		FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);

	// scale output:
	double reciprocalSize = 1.0 / input.size();
	for (auto &c : output){
		c *= reciprocalSize;
	}

	return output;
}

// AnalyticSignalV() : Analytic signal of vector of Complex doubles
// (Note: This function is referred to as "hilbert()" in Matlab / Octave, but it is not exactly a hilbert transform. 
// The hilbert Transform is placed in the imaginary part, and the original input is in the real part.)
// See Footnote* below for more information on algorithm ...

std::vector<std::complex<double>>
AnalyticSignalV(const std::vector<std::complex<double>>& input) {

	std::vector<std::complex<double>> U = fftV(input);

	size_t N = input.size();
	size_t halfN = N / 2;
	
	// Note: U[0], U[halfN] unchanged:
	for (size_t n = 1; n < N; ++n) {
		if (n > halfN)
			U[n] = 0;
		if (n < halfN)
			U[n] *= 2.0;
	}
	
	std::vector<std::complex<double>> output = ifftV(U);
	return output;
}

// makeMinPhase() : transform linear-phase FIR filter coefficients into minimum-phase (in-place version)
template<typename FloatType>
void makeMinPhase(FloatType* pFIRcoeffs, size_t length)
{
	size_t fftLength = pow(2, 2.0 + ceil(log2(length))); // use FFT 4x larger than (length rounded-up to power-of-2)

	std::vector <std::complex<double>> complexInput;
	std::vector <std::complex<double>> complexOutput;

	// Pad zeros on either side of FIR:
	
	size_t frontPaddingLength = (fftLength - length) / 2;
	size_t backPaddingLength = fftLength - frontPaddingLength - length;
	
	for (size_t n = 0; n < frontPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	for (size_t n = 0; n < length; ++n) {
		complexInput.push_back({ pFIRcoeffs[n], 0 });
	}

	for (size_t n = 0; n < backPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	/*
	// pad with trailing zeros
	for (int n = 0; n < fftLength; ++n) {
		if (n<length)
			complexInput.push_back({ pFIRcoeffs[n], 0 });
		else
			complexInput.push_back({ 0, 0 }); // pad remainder with zeros
	}
	*/

	assert(complexInput.size() == fftLength); // make sure padding worked properly.

	// Formula is as follows:

	// take the reversed array of
	// the real parts of
	// the ifft of
	// e to the power of
	// the Analytic Signal of
	// the real parts of 
	// the log of
	// the dynamic-ranged limited version of
	// the fft of 
	// the original filter
			
	complexOutput = realV(ifftV(expV(AnalyticSignalV(realV(logV(limitDynRangeV(fftV(complexInput),-190)))))));
	std::reverse(complexOutput.begin(), complexOutput.end());

	// write all the real parts back to coeff array:
	int n = 0;
	for (auto &c : complexOutput) {
		if (n < length)
			pFIRcoeffs[n] = c.real();	
		else
			break;
		++n;
	}
}

// makeMinPhase2() : take linear-phase FIR filter coefficients, and return a new vector of minimum-phase coefficients
template<typename FloatType>
std::vector<FloatType> makeMinPhase2(const FloatType* pFIRcoeffs, size_t length)
{
	size_t fftLength = pow(2, 2.0 + ceil(log2(length))); // use FFT 4x larger than (length rounded-up to power-of-2)

	std::vector <std::complex<double>> complexInput;
	std::vector <std::complex<double>> complexOutput;

	// Pad zeros on either side of FIR:

	size_t frontPaddingLength = (fftLength - length) / 2;
	size_t backPaddingLength = fftLength - frontPaddingLength - length;

	for (size_t n = 0; n < frontPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	for (size_t n = 0; n < length; ++n) {
		complexInput.push_back({ pFIRcoeffs[n], 0 });
	}

	for (size_t n = 0; n < backPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	assert(complexInput.size() == fftLength); // make sure padding worked properly.

	// Formula is as follows:

	// take the reversed array of
	// the real parts of
	// the ifft of
	// e to the power of
	// the Analytic Signal of
	// the real parts of 
	// the log of
	// the dynamic-ranged limited version of
	// the fft of 
	// the original filter

	complexOutput = realV(ifftV(expV(AnalyticSignalV(realV(logV(limitDynRangeV(fftV(complexInput), -190)))))));
	std::reverse(complexOutput.begin(), complexOutput.end());

	std::vector<FloatType> minPhaseCoeffs;
	minPhaseCoeffs.reserve(complexOutput.size());
	for (auto & c : complexOutput) {
		minPhaseCoeffs.push_back(c.real());
	}

	return minPhaseCoeffs;
}

///////////////////////////////////////////////////////////////////////
// utility functions:

// dumpKaiserWindow() - utility function for displaying Kaiser Window:
void dumpKaiserWindow(size_t length, double Beta) {
    std::vector<double> f(length, 1);
    applyKaiserWindow<double>(f.data(), length, Beta);
    for (int i = 0; i < length; ++i) {
        std::cout << i << ": " << f[i] << std::endl;
    }

    std::vector<double> g(length, 1);
    applyKaiserWindow<double>(g.data(), length, Beta);
    for (int i = 0; i < length; ++i) {
        std::cout << i << ": " << g[i] << std::endl;
    }
}

// asserts that the two Kaiser Window formulas agree with each other (within a specified tolerance)
void assertKaiserWindow(size_t length, double Beta) {

	const double tolerance = 0.001;
	const double upper = 1.0 + tolerance;
	const double lower = 1.0 - tolerance;

	std::vector<double> f(length, 1);
	applyKaiserWindow2<double>(f.data(), length, Beta);

	std::vector<double> g(length, 1);
	applyKaiserWindow<double>(g.data(), length, Beta);

	for (int i = 0; i < length; ++i) {
		double ratio = f[i] / g[i];
		assert(ratio < upper && ratio > lower);
	}
}

// dumpFilter() - utility function for displaying filter coefficients:
template<typename FloatType> void dumpFilter(const FloatType* Filter, int Length) {
	for (int i = 0; i < Length; ++i) {
		std::cout << Filter[i] << std::endl;
	}
}

void testMinPhase() {
	std::cout << std::setprecision(15);
	const size_t MedFilterSize = 511;
	double MedFilterTaps[511];
	makeLPF<double>(MedFilterTaps, MedFilterSize, 21819, 96000);
	applyKaiserWindow2<double>(MedFilterTaps, MedFilterSize, calcKaiserBeta(195));
	dumpFilter(MedFilterTaps, MedFilterSize);
	makeMinPhase(MedFilterTaps, MedFilterSize);
	dumpFilter(MedFilterTaps, MedFilterSize);
}

void dumpComplexVector(const std::vector<std::complex<double>>& v)
{
	for (auto &c : v) {
		std::cout << c.real() << "+" << c.imag() << "i" << std::endl;
	}
}

template<typename FloatType>
void dumpFFT(FloatType* data, size_t length)
{
	size_t pow2length = pow(2, 1.0 + floor(log2(length)));

	std::vector <std::complex<double>> complexInput;
	std::vector <std::complex<double>> complexOutput;

	for (int n = 0; n < pow2length; ++n) {
		if (n<length)
			complexInput.push_back({ data[n], 0 });
		else
			complexInput.emplace_back(0, 0); // pad remainder with zeros (to-do: does it mattter where the zeros are put ?)
	}
	
	complexOutput = fftV(complexInput);

	std::setprecision(17);
	std::cout << "real,imag,mag,phase" << std::endl;
	for (auto &c : complexOutput) {
		std::cout << c.real() << "," << c.imag() << "," << std::abs(c) << "," << arg(c) << std::endl;
	}
}

void testSinAccuracy() {
	
	const int numSteps = 10000000;
	const double inc = M_PI / numSteps;
	double t = M_PI / -2.0;
	double maxError = 0.0;
	double worstT = 0.0;

	for (int i = 0; i < numSteps; ++i ) {
		// calc relative error of
		// |(sin 2t - 2 * sint * cost) / sin 2t|
		// (double-angle identity)
		
		double e = std::abs((std::sin(2.0 * t) - 2.0 * std::sin(t) * std::cos(t)) / std::sin(2.0 * t));
		//double e = std::abs((sin(2.0 * t) - 2.0 * sin(t) * cos(t)) / sin(2.0 * t));
		if (e > maxError) {
			worstT = t;
			maxError = e;
		}
		t += inc;
	}
	std::cout << "maxError: " << std::setprecision(33) << maxError << std::endl;
	std::cout << "worstT: " << worstT << std::endl;
}

// *Marple, S. L. "Computing the Discrete-Time Analytic Signal via FFT." IEEE Transactions on Signal Processing. Vol. 47, 1999, pp. 2600�2603

#endif // FIRFFILTER_H_