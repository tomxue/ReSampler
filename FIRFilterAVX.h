/*
* Copyright (C) 2016 - 2018 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef FIRFFILTER_AVX_H_
#define FIRFFILTER_AVX_H_

// FIRFilterAVX.h : AVX-specific code for FIR filtering

#include <immintrin.h>
#include <typeinfo>
#include <cstdint>
#include <cassert>
#include "alignedmalloc.h"

#ifdef USE_AVX
#define AVX_ALIGNMENT_SIZE 32
static inline float sum8floats(__m256 x);
static inline double sum4doubles(__m256d x);

//#define USE_FMA 1

template <typename FloatType>
class FIRFilter {

public:
	// constructor:
	FIRFilter(const FloatType* taps, size_t size) :
		size(size), sizeRounded8((size >> 3) << 3), sizeRounded4((size >> 2) << 2), currentIndex(size-1), lastPut(0),
		signal(nullptr), kernel0(nullptr), kernel1(nullptr), kernel2(nullptr), kernel3(nullptr), kernel4(nullptr), kernel5(nullptr), kernel6(nullptr), kernel7(nullptr)
	{
		// allocate buffers:
		allocateBuffers();
		assertAlignment();

		// initialize filter kernel and signal buffers
		for (unsigned int i = 0; i < size; ++i) {
			kernel0[i] = taps[i];
			signal[i] = 0.0;
			signal[i + size] = 0.0;
		}

		// Populate additional kernel Phases:
		memcpy(1 + kernel1, kernel0, (size - 1) * sizeof(FloatType));
		kernel1[0] = kernel0[size - 1];
		memcpy(1 + kernel2, kernel1, (size - 1) * sizeof(FloatType));
		kernel2[0] = kernel1[size - 1];
		memcpy(1 + kernel3, kernel2, (size - 1) * sizeof(FloatType));
		kernel3[0] = kernel2[size - 1];
		//
		memcpy(1 + kernel4, kernel3, (size - 1) * sizeof(FloatType));
		kernel4[0] = kernel3[size - 1];
		memcpy(1 + kernel5, kernel4, (size - 1) * sizeof(FloatType));
		kernel5[0] = kernel4[size - 1];
		memcpy(1 + kernel6, kernel5, (size - 1) * sizeof(FloatType));
		kernel6[0] = kernel5[size - 1];
		memcpy(1 + kernel7, kernel6, (size - 1) * sizeof(FloatType));
		kernel7[0] = kernel6[size - 1];

	}

	// deconstructor:
	~FIRFilter() {
		freeBuffers();
	}

	// copy constructor: 
	FIRFilter(const FIRFilter& other) : size(other.size), sizeRounded8(other.sizeRounded8), sizeRounded4(other.sizeRounded4), currentIndex(other.currentIndex), lastPut(other.lastPut)
	{
		allocateBuffers();
		assertAlignment();
		copyBuffers(other);
	}

	// move constructor:
	FIRFilter(FIRFilter&& other) :
		size(other.size), sizeRounded8(other.sizeRounded8), sizeRounded4(other.sizeRounded4), currentIndex(other.currentIndex), lastPut(other.lastPut),
		signal(other.signal), kernel0(other.kernel0), kernel1(other.kernel1), kernel2(other.kernel2), kernel3(other.kernel3),
		kernel4(other.kernel4), kernel5(other.kernel5), kernel6(other.kernel6), kernel7(other.kernel7)
	{
		assertAlignment();
		other.signal = nullptr;
		other.kernel0 = nullptr;
		other.kernel1 = nullptr;
		other.kernel2 = nullptr;
		other.kernel3 = nullptr;
		other.kernel4 = nullptr;
		other.kernel5 = nullptr;
		other.kernel6 = nullptr;
		other.kernel7 = nullptr;
	}

	// copy assignment:
	FIRFilter& operator= (const FIRFilter& other)
	{
		size = other.size;
		sizeRounded8 = other.sizeRounded8;
		sizeRounded4 = other.sizeRounded4;
		currentIndex = other.currentIndex;
		lastPut = other.lastPut;
		freeBuffers();
		allocateBuffers();
		assertAlignment();
		copyBuffers(other);
		return *this;
	}

	// move assignment:
	FIRFilter& operator= (FIRFilter&& other)
	{
		if (this != &other) // prevent self-assignment
		{
			size = other.size;
			sizeRounded8 = other.sizeRounded8;
			sizeRounded4 = other.sizeRounded4;
			currentIndex = other.currentIndex;
			lastPut = other.lastPut;

			freeBuffers();

			signal = other.signal;
			kernel0 = other.kernel0;
			kernel1 = other.kernel1;
			kernel2 = other.kernel2;
			kernel3 = other.kernel3;
			kernel4 = other.kernel4;
			kernel5 = other.kernel5;
			kernel6 = other.kernel6;
			kernel7 = other.kernel7;

			assertAlignment();

			other.signal = nullptr;
			other.kernel0 = nullptr;
			other.kernel1 = nullptr;
			other.kernel2 = nullptr;
			other.kernel3 = nullptr;
			other.kernel4 = nullptr;
			other.kernel5 = nullptr;
			other.kernel6 = nullptr;
			other.kernel7 = nullptr;
		}
		return *this;
	}
	
	void reset() {
		// reset indexes:
		currentIndex = size - 1;
		lastPut = 0;

		// clear signal buffer
		for (unsigned int i = 0; i < size; ++i) {
			signal[i] = 0.0;
			signal[i + size] = 0.0;
		}

	}

	void put(FloatType value) { // Put signal in reverse order.
		signal[currentIndex] = value;
		lastPut = currentIndex;
		if (currentIndex == 0) {
			currentIndex = size - 1; // Wrap
			memcpy(signal + size, signal, size*sizeof(FloatType)); // copy history to upper half of buffer
		}
		else
			--currentIndex;
	}

	void putZero() {
		signal[currentIndex] = 0.0;
		if (currentIndex == 0) {
			currentIndex = size - 1; // Wrap
			memcpy(signal + size, signal, size*sizeof(FloatType)); // copy history to upper half of buffer
		}
		else
			--currentIndex;
	}

	FloatType get() {

		// AVX implementation: This only works with floats - doubles need specialization ...

		FloatType output = 0.0;
		FloatType* kernel = kernel0;
		int index = (currentIndex >> 3) << 3; // make multiple-of-eight
		int phase = currentIndex & 7;
		
		// Part1 : Head
		// select proper Kernel phase and calculate first Block of 8:
		switch (phase) {

		case 0:
			kernel = kernel0;
			output =
				kernel[0] * signal[index] + kernel[1] * signal[index + 1] + kernel[2] * signal[index + 2] + kernel[3] * signal[index + 3] +
				kernel[4] * signal[index + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 1:
			kernel = kernel1;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + 1] + kernel[2] * signal[index + 2] + kernel[3] * signal[index + 3] +
				kernel[4] * signal[index + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 2:
			kernel = kernel2;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + 2] + kernel[3] * signal[index + 3] +
				kernel[4] * signal[index + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 3: 
			kernel = kernel3;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + size + 2] + kernel[3] * signal[index + 3] +
				kernel[4] * signal[index + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 4:
			kernel = kernel4;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + size + 2] + kernel[3] * signal[index + size + 3] +
				kernel[4] * signal[index + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 5:
			kernel = kernel5;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + size + 2] + kernel[3] * signal[index + size + 3] +
				kernel[4] * signal[index + size + 4] + kernel[5] * signal[index + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 6:
			kernel = kernel6;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + size + 2] + kernel[3] * signal[index + size + 3] +
				kernel[4] * signal[index + size + 4] + kernel[5] * signal[index + size + 5] + kernel[6] * signal[index + 6] + kernel[7] * signal[index + 7];
			break;

		case 7:
			kernel = kernel7;
			output =
				kernel[0] * signal[index + size] + kernel[1] * signal[index + size + 1] + kernel[2] * signal[index + size + 2] + kernel[3] * signal[index + size + 3] +
				kernel[4] * signal[index + size + 4] + kernel[5] * signal[index + size + 5] + kernel[6] * signal[index + size + 6] + kernel[7] * signal[index + 7];
			break;
		}
		index += 8;

		// Part 2: Body
		alignas(AVX_ALIGNMENT_SIZE) __m256 s;	// AVX Vector Registers for calculation
		alignas(AVX_ALIGNMENT_SIZE) __m256 k;
		alignas(AVX_ALIGNMENT_SIZE) __m256 product;
		alignas(AVX_ALIGNMENT_SIZE) __m256 accumulator = _mm256_setzero_ps();

		for (int i = 8; i < sizeRounded8; i += 8) {
			s = _mm256_load_ps(signal + index);
			k = _mm256_load_ps(kernel + i);
#ifdef USE_FMA
			accumulator = _mm256_fmadd_ps(signal, kernel, accumulator);
#else
			product = _mm256_mul_ps(s, k);
			accumulator = _mm256_add_ps(product, accumulator);
#endif

			index += 8;
		}

		output += sum8floats(accumulator);

		// Part 3: Tail
		for (int j = sizeRounded8; j < size; ++j) {
			output += signal[index] * kernel[j];
			++index;
		}

		return output;
	}

	FloatType lazyGet(int L) {	// Skips stuffed-zeros introduced by interpolation, by only calculating every Lth sample from lastPut
		FloatType output = 0.0;
		int Offset = lastPut - currentIndex;
		if (Offset < 0) { // Wrap condition
			Offset += size;
		}
	
		for (int i = Offset; i < size; i+=L) {
			output += signal[i+ currentIndex] * kernel0[i];
		}
		return output;
	}

private:
	size_t size;
	size_t sizeRounded8; // size rounded down to nearest multiple of 8
	size_t sizeRounded4; // size rounded down to nearest multiple of 4
	FloatType* signal; // Double-length signal buffer, to facilitate fast emulation of a circular buffer
	int currentIndex;
	int lastPut;

	// Polyphase Filter Kernel table:
	FloatType* kernel0;
	FloatType* kernel1;
	FloatType* kernel2;
	FloatType* kernel3;
	FloatType* kernel4;
	FloatType* kernel5;
	FloatType* kernel6;
	FloatType* kernel7;

	void allocateBuffers()
	{
		signal = static_cast<FloatType*>(aligned_malloc(2 * size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel0 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel1 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel2 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel3 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel4 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel5 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel6 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
		kernel7 = static_cast<FloatType*>(aligned_malloc(size * sizeof(FloatType), AVX_ALIGNMENT_SIZE));
	}

	void copyBuffers(const FIRFilter& other)
	{
		memcpy(signal, other.signal, 2 * size * sizeof(FloatType));
		memcpy(kernel0, other.kernel0, size * sizeof(FloatType));
		memcpy(kernel1, other.kernel1, size * sizeof(FloatType));
		memcpy(kernel2, other.kernel2, size * sizeof(FloatType));
		memcpy(kernel3, other.kernel3, size * sizeof(FloatType));
		memcpy(kernel4, other.kernel4, size * sizeof(FloatType));
		memcpy(kernel5, other.kernel5, size * sizeof(FloatType));
		memcpy(kernel6, other.kernel6, size * sizeof(FloatType));
		memcpy(kernel7, other.kernel7, size * sizeof(FloatType));
	}

	void freeBuffers()
	{
		aligned_free(signal);
		aligned_free(kernel0);
		aligned_free(kernel1);
		aligned_free(kernel2);
		aligned_free(kernel3);
		aligned_free(kernel4);
		aligned_free(kernel5);
		aligned_free(kernel6);
		aligned_free(kernel7);
	}

	// assertAlignment() : asserts that all private data buffers are aligned on expected boundaries
	void assertAlignment()
	{
		const std::uintptr_t alignment = AVX_ALIGNMENT_SIZE;
		
		assert(reinterpret_cast<std::uintptr_t>(signal) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel0) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel1) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel2) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel3) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel4) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel5) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel6) % alignment == 0);
		assert(reinterpret_cast<std::uintptr_t>(kernel7) % alignment == 0);

	}
};

// ================================= 
// AVX specializations for doubles :
// =================================

template <>
double FIRFilter<double>::get() {

	// AVX implementation: This only works with doubles !
	// Processes four doubles at a time.

	double output = 0.0;
	double* Kernel;
	int index = (currentIndex >> 2) << 2; // make multiple-of-four
	int phase = currentIndex & 3;

	// Part1 : Head
	// select proper Kernel phase and calculate first Block of 4:
	switch (phase) {

	case 0:
		Kernel = kernel0;
		output = Kernel[0] * signal[index] + Kernel[1] * signal[index + 1] + Kernel[2] * signal[index + 2] + Kernel[3] * signal[index + 3];
		break;

	case 1:
		Kernel = kernel1;
		output = Kernel[0] * signal[index + size] + Kernel[1] * signal[index + 1] + Kernel[2] * signal[index + 2] + Kernel[3] * signal[index + 3];
		break;

	case 2:
		Kernel = kernel2;
		output = Kernel[0] * signal[index + size] + Kernel[1] * signal[index + size + 1] + Kernel[2] * signal[index + 2] + Kernel[3] * signal[index + 3];
		break;

	case 3:
		Kernel = kernel3;
		output = Kernel[0] * signal[index + size] + Kernel[1] * signal[index + size + 1] + Kernel[2] * signal[index + size + 2] + Kernel[3] * signal[index + 3];
		break;

	}
	
	index += 4;

	// Part 2: Body
	alignas(AVX_ALIGNMENT_SIZE) __m256d s;	// AVX Vector Registers for calculation
	alignas(AVX_ALIGNMENT_SIZE) __m256d k;
	alignas(AVX_ALIGNMENT_SIZE) __m256d product;
	alignas(AVX_ALIGNMENT_SIZE) __m256d accumulator = _mm256_setzero_pd();

	for (int i = 4; i < sizeRounded4; i += 4) {
		s = _mm256_load_pd(signal + index);
		k = _mm256_load_pd(Kernel + i);

#ifdef USE_FMA
		accumulator = _mm256_fmadd_pd(signal, kernel, accumulator);
#else
		product = _mm256_mul_pd(s, k);
		accumulator = _mm256_add_pd(product, accumulator);
#endif
		index += 4;
	}

	output += sum4doubles(accumulator);

	// Part 3: Tail
	for (int j = sizeRounded4; j < size; ++j) {
		output += signal[index] * Kernel[j];
		++index;
	}
	return output;
}

// Horizontal add function (sums 8 floats into single float) http://stackoverflow.com/questions/23189488/horizontal-sum-of-32-bit-floats-in-256-bit-avx-vector
static inline float sum8floats(__m256 x) {
	const __m128 x128 = _mm_add_ps(
		_mm256_extractf128_ps(x, 1),
		_mm256_castps256_ps128(x));																// ( x3+x7, x2+x6, x1+x5, x0+x4 )
	const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));								// ( -, -, x1+x3+x5+x7, x0+x2+x4+x6 )
	const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));							// ( -, -, -, x0+x1+x2+x3+x4+x5+x6+x7 )
	return _mm_cvtss_f32(x32);
}

// Horizontal add function (sums 4 doubles into single double)
static inline double sum4doubles(__m256d x) {
	const __m128d x128 = _mm_add_pd(
		_mm256_extractf128_pd(x, 1),
		_mm256_castpd256_pd128(x));
	const __m128d x64 = _mm_add_pd(_mm_permute_pd(x128, 1), x128);
	return _mm_cvtsd_f64(x64);
}

double testSum4Doubles(double a, double b, double c, double d) {
	const __m256d v = _mm256_set_pd(a, b, c, d);
	return sum4doubles(v);
}

#endif // USE_AVX
#endif // FIRFFILTER_AVX_H_