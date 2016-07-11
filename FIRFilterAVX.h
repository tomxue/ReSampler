#ifndef FIRFFILTER_AVX_H_
#define FIRFFILTER_AVX_H_

//#define USE_FMA 1

// FIRFilterAVX.h : AVX-specific code for FIR filtering

#include <immintrin.h>
#include <typeinfo>

#define FILTERSIZE_HUGE 32767
#define FILTERSIZE_MEDIUM 511

template <typename FloatType, unsigned int size>
class FIRFilter {
public:

	FIRFilter(FloatType* taps) :
		CurrentIndex(size-1), LastPut(0)
	{
		for (unsigned int i = 0; i < size; ++i) {
			Kernel0[i] = taps[i];
			Signal[i] = 0.0;
			Signal[i + size] = 0.0;
		}

#ifdef USE_AVX
		// Populate remaining kernel Phases:
		memcpy(1 + Kernel1, Kernel0, (size - 1) * sizeof(FloatType));
		Kernel1[0] = Kernel0[size - 1];
		memcpy(1 + Kernel2, Kernel1, (size - 1) * sizeof(FloatType));
		Kernel2[0] = Kernel1[size - 1];
		memcpy(1 + Kernel3, Kernel2, (size - 1) * sizeof(FloatType));
		Kernel3[0] = Kernel2[size - 1];
		//
		memcpy(1 + Kernel4, Kernel3, (size - 1) * sizeof(FloatType));
		Kernel4[0] = Kernel3[size - 1];
		memcpy(1 + Kernel5, Kernel4, (size - 1) * sizeof(FloatType));
		Kernel5[0] = Kernel4[size - 1];
		memcpy(1 + Kernel6, Kernel5, (size - 1) * sizeof(FloatType));
		Kernel6[0] = Kernel5[size - 1];
		memcpy(1 + Kernel7, Kernel6, (size - 1) * sizeof(FloatType));
		Kernel7[0] = Kernel6[size - 1];

#endif

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

#ifndef USE_AVX
		FloatType output = 0.0;
		int index = CurrentIndex;
		for (int i = 0; i < size; ++i) {
			output += Signal[index] * Kernel0[i];
			index++;
		}
		return output;
#else
		// AVX implementation: This only works with floats !

		FloatType output = 0.0;
		FloatType* Kernel = Kernel0;
		int Index = (CurrentIndex >> 3) << 3; // make multiple-of-eight
		int Phase = CurrentIndex & 7;
		
		// Part1 : Head
		// select proper Kernel phase and calculate first Block of 8:
		switch (Phase) {

		case 0:
			Kernel = Kernel0;
			output =
				Kernel[0] * Signal[Index] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3] +
				Kernel[4] * Signal[Index + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 1:
			Kernel = Kernel1;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3] +
				Kernel[4] * Signal[Index + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 2:
			Kernel = Kernel2;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3] +
				Kernel[4] * Signal[Index + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 3: 
			Kernel = Kernel3;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + 3] +
				Kernel[4] * Signal[Index + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 4:
			Kernel = Kernel4;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + size + 3] +
				Kernel[4] * Signal[Index + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 5:
			Kernel = Kernel5;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + size + 3] +
				Kernel[4] * Signal[Index + size + 4] + Kernel[5] * Signal[Index + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 6:
			Kernel = Kernel6;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + size + 3] +
				Kernel[4] * Signal[Index + size + 4] + Kernel[5] * Signal[Index + size + 5] + Kernel[6] * Signal[Index + 6] + Kernel[7] * Signal[Index + 7];
			break;

		case 7:
			Kernel = Kernel7;
			output =
				Kernel[0] * Signal[Index + size] + Kernel[1] * Signal[Index + size + 1] + Kernel[2] * Signal[Index + size + 2] + Kernel[3] * Signal[Index + size + 3] +
				Kernel[4] * Signal[Index + size + 4] + Kernel[5] * Signal[Index + size + 5] + Kernel[6] * Signal[Index + size + 6] + Kernel[7] * Signal[Index + 7];
			break;
		}
		Index += 8;

		// Part 2: Body
		alignas(32) __m256 signal;	// AVX Vector Registers for calculation
		alignas(32) __m256 kernel;
		alignas(32) __m256 product;
		alignas(32) __m256 accumulator = _mm256_setzero_ps();

		for (int i = 8; i < (size >> 3) << 3; i += 8) {
			signal = _mm256_load_ps(Signal + Index);
			kernel = _mm256_load_ps(Kernel + i);
#ifdef USE_FMA
			accumulator = _mm256_fmadd_ps(signal, kernel, accumulator);
#else
			product = _mm256_mul_ps(signal, kernel);
			accumulator = _mm256_add_ps(product, accumulator);
#endif

			Index += 8;
		}	

		output += 
			accumulator.m256_f32[0] +
			accumulator.m256_f32[1] +
			accumulator.m256_f32[2] +
			accumulator.m256_f32[3] +
			accumulator.m256_f32[4] +
			accumulator.m256_f32[5] +
			accumulator.m256_f32[6] +
			accumulator.m256_f32[7];

		// Part 3: Tail
		for (int j = (size >> 3) << 3; j < size; ++j) {
			output += Signal[Index] * Kernel[j];
			++Index;
		}

		return output;

#endif // !USE_AVX
	}

	FloatType LazyGet(int L) {	// Skips stuffed-zeros introduced by interpolation, by only calculating every Lth sample from LastPut
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
	alignas(32) FloatType Kernel0[size];

#ifdef USE_AVX
	// Polyphase Filter Kernel table:
	alignas(32) FloatType Kernel1[size];
	alignas(32) FloatType Kernel2[size];
	alignas(32) FloatType Kernel3[size];
	alignas(32) FloatType Kernel4[size];
	alignas(32) FloatType Kernel5[size];
	alignas(32) FloatType Kernel6[size];
	alignas(32) FloatType Kernel7[size];

#endif
	alignas(32) FloatType Signal[size*2];	// Double-length signal buffer, to facilitate fast emulation of a circular buffer
	int CurrentIndex;
	int LastPut;
};

#ifdef USE_AVX

// ================================= 
// AVX specializations for doubles :
// =================================

double FIRFilter<double, FILTERSIZE_MEDIUM>::get() {

	// AVX implementation: This only works with doubles !
	// Processes four doubles at a time.

	double output = 0.0;
	double* Kernel;
	int Index = (CurrentIndex >> 2) << 2; // make multiple-of-four
	int Phase = CurrentIndex & 3;

	// Part1 : Head
	// select proper Kernel phase and calculate first Block of 4:
	switch (Phase) {

	case 0:
		Kernel = Kernel0;
		output = Kernel[0] * Signal[Index] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 1:
		Kernel = Kernel1;
		output = Kernel[0] * Signal[Index + FILTERSIZE_MEDIUM] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 2:
		Kernel = Kernel2;
		output = Kernel[0] * Signal[Index + FILTERSIZE_MEDIUM] + Kernel[1] * Signal[Index + FILTERSIZE_MEDIUM + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 3:
		Kernel = Kernel3;
		output = Kernel[0] * Signal[Index + FILTERSIZE_MEDIUM] + Kernel[1] * Signal[Index + FILTERSIZE_MEDIUM + 1] + Kernel[2] * Signal[Index + FILTERSIZE_MEDIUM + 2] + Kernel[3] * Signal[Index + 3];
		break;

	}
	Index += 4;

	// Part 2: Body
	alignas(32) __m256d signal;	// AVX Vector Registers for calculation
	alignas(32) __m256d kernel;
	alignas(32) __m256d product;
	alignas(32) __m256d accumulator = _mm256_setzero_pd();

	for (int i = 4; i < (FILTERSIZE_MEDIUM >> 2) << 2; i += 4) {
		signal = _mm256_load_pd(Signal + Index);
		kernel = _mm256_load_pd(Kernel + i);
		product = _mm256_mul_pd(signal, kernel);
		accumulator = _mm256_add_pd(product, accumulator);
		Index += 4;
	}

	output +=
		accumulator.m256d_f64[0] +
		accumulator.m256d_f64[1] +
		accumulator.m256d_f64[2] +
		accumulator.m256d_f64[3];

	// Part 3: Tail
	for (int j = (FILTERSIZE_MEDIUM >> 2) << 2; j < FILTERSIZE_MEDIUM; ++j) {
		output += Signal[Index] * Kernel[j];
		++Index;
	}

	return output;
}

double FIRFilter<double, FILTERSIZE_HUGE>::get() {

	// AVX implementation: This only works with doubles !
	// Processes four doubles at a time.

	double output = 0.0;
	double* Kernel;
	int Index = (CurrentIndex >> 2) << 2; // make multiple-of-four
	int Phase = CurrentIndex & 3;

	// Part1 : Head
	// select proper Kernel phase and calculate first Block of 4:
	switch (Phase) {

	case 0:
		Kernel = Kernel0;
		output = Kernel[0] * Signal[Index] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 1:
		Kernel = Kernel1;
		output = Kernel[0] * Signal[Index + FILTERSIZE_HUGE] + Kernel[1] * Signal[Index + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 2:
		Kernel = Kernel2;
		output = Kernel[0] * Signal[Index + FILTERSIZE_HUGE] + Kernel[1] * Signal[Index + FILTERSIZE_HUGE + 1] + Kernel[2] * Signal[Index + 2] + Kernel[3] * Signal[Index + 3];
		break;

	case 3:
		Kernel = Kernel3;
		output = Kernel[0] * Signal[Index + FILTERSIZE_HUGE] + Kernel[1] * Signal[Index + FILTERSIZE_HUGE + 1] + Kernel[2] * Signal[Index + FILTERSIZE_HUGE + 2] + Kernel[3] * Signal[Index + 3];
		break;

	}
	Index += 4;

	// Part 2: Body
	alignas(32) __m256d signal;	// AVX Vector Registers for calculation
	alignas(32) __m256d kernel;
	alignas(32) __m256d product;
	alignas(32) __m256d accumulator = _mm256_setzero_pd();

	for (int i = 4; i < (FILTERSIZE_HUGE >> 2) << 2; i += 4) {
		signal = _mm256_load_pd(Signal + Index);
		kernel = _mm256_load_pd(Kernel + i);
		product = _mm256_mul_pd(signal, kernel);
		accumulator = _mm256_add_pd(product, accumulator);
		Index += 4;
	}

	output +=
		accumulator.m256d_f64[0] +
		accumulator.m256d_f64[1] +
		accumulator.m256d_f64[2] +
		accumulator.m256d_f64[3];

	// Part 3: Tail
	for (int j = (FILTERSIZE_HUGE >> 2) << 2; j < FILTERSIZE_HUGE; ++j) {
		output += Signal[Index] * Kernel[j];
		++Index;
	}

	return output;
}
#endif // USE_AVX

#endif // FIRFFILTER_AVX_H_