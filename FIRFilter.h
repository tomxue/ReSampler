#ifndef FIRFFILTER_H_
#define FIRFFILTER_H_

#include <typeinfo>

#define FILTERSIZE_HUGE 32767
#define FILTERSIZE_MEDIUM 511


#define USE_SIMD 1 // 2016/04/01: Needs specializations (SIMD code won't work for double precision)

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

#ifdef USE_SIMD
		// Populate remaining kernel Phases:
		memcpy(1 + Kernel1, Kernel0, (size - 1)*sizeof(FloatType));
		Kernel1[0] = Kernel0[size - 1];
		memcpy(1 + Kernel2, Kernel1, (size - 1)*sizeof(FloatType));
		Kernel2[0] = Kernel1[size - 1];
		memcpy(1 + Kernel3, Kernel2, (size - 1)*sizeof(FloatType));
		Kernel3[0] = Kernel2[size - 1];
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

#ifndef USE_SIMD
		FloatType output = 0.0;
		int index = CurrentIndex;
		for (int i = 0; i < size; ++i) {
			output += Signal[index] * Kernel0[i];
			index++;
		}
		return output;
#else
		// SIMD implementation: This only works with floats !

		FloatType output = 0.0;
		FloatType* Kernel;
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
		alignas(16) __m128 signal;	// SIMD Vector Registers for calculation
		alignas(16) __m128 kernel;
		alignas(16) __m128 product;
		alignas(16) __m128 accumulator = _mm_setzero_ps();

		for (int i = 4; i < (size >> 2) << 2; i += 4) {
			signal = _mm_load_ps(Signal + Index);
			kernel = _mm_load_ps(Kernel + i);
			product = _mm_mul_ps(signal, kernel);
			accumulator = _mm_add_ps(product, accumulator);
			Index += 4;
		}	

		output += accumulator.m128_f32[0] +
			accumulator.m128_f32[1] +
			accumulator.m128_f32[2] +
			accumulator.m128_f32[3];

		// Part 3: Tail
		for (int j = (size >> 2) << 2; j < size; ++j) {
			output += Signal[Index] * Kernel[j];
			++Index;
		}

		return output;

#endif // !USE_SIMD
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
	alignas(16) FloatType Kernel0[size];

#ifdef USE_SIMD
	// Polyphase Filter Kernel table:
	alignas(16) FloatType Kernel1[size];
	alignas(16) FloatType Kernel2[size];
	alignas(16) FloatType Kernel3[size];
#endif

	alignas(16) FloatType Signal[size*2];	// Double-length signal buffer, to facilitate fast emulation of a circular buffer
	int CurrentIndex;
	int LastPut;
};

#ifdef USE_SIMD

// Specializations for Doubles:

double FIRFilter<double, FILTERSIZE_MEDIUM>::get() {
	double output = 0.0;
	int index = CurrentIndex;
	for (int i = 0; i < FILTERSIZE_MEDIUM; ++i) {
		output += Signal[index] * Kernel0[i];
		index++;
	}
	return output;
}

double FIRFilter<double, FILTERSIZE_HUGE>::get() {
	double output = 0.0;
	int index = CurrentIndex;
	for (int i = 0; i < FILTERSIZE_HUGE; ++i) {
		output += Signal[index] * Kernel0[i];
		index++;
	}
	return output;
}
#endif // USE_SIMD

#endif // FIRFFILTER_H_
