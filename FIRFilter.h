#ifndef FIRFFILTER_H_
#define FIRFFILTER_H_

//#define USE_SIMD 1 // 13/03/2016; Still NQR !

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
		FloatType output = 0.0;
		int index = CurrentIndex;
		for (int i = 0; i < size; ++i) {
			output += Signal[index] * Kernel0[i];
			index++;
		}
		return output;
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
	alignas(16) FloatType Signal[size*2];	// Double-length signal buffer, to facilitate fast emulation of a circular buffer
	int CurrentIndex;
	int LastPut;
};

#endif