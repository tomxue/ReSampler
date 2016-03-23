#ifndef FIRFilter_H_
#define FIRFFILTER_H_

// #define USE_SIMD 1 // 13/03/2016; Still NQR !

template <typename FloatType, unsigned int size>
class FIRFilter {
public:

	FIRFilter(FloatType* taps) :
		m_CurrentIndex(0) 
	{
		for (unsigned int i = 0; i < size; ++i) {
			m_Taps[i] = taps[i];
			m_Signal[i] = 0.0;
		}
	}

	void put(FloatType value) {
		m_Signal[m_CurrentIndex++] = value;
		if (m_CurrentIndex == size)
			m_CurrentIndex = 0;
	}

	void putZero() {
		m_Signal[m_CurrentIndex++] = 0.0;
		if (m_CurrentIndex == size)
			m_CurrentIndex = 0;
	}

	FloatType get() {

		FloatType output = 0.0;

#ifdef USE_SIMD

		// Reverse and align signal:
		int index = m_CurrentIndex;
		for (int i = 0; i < size; ++i) {
			m_ReversedSignal[i] = m_Signal[index];
			index = (index == 0) ? size - 1 : index - 1;
		}

		int MultOf4Size = (size >> 2) << 2; // nearest multiple-of-4 below size

		alignas(16) __m128 signal;
		alignas(16) __m128 kernel;
		alignas(16) __m128 product;
		alignas(16) __m128 accumulator = _mm_setzero_ps();

		for (int i = 0; i < MultOf4Size; i += 4) {
			signal = _mm_load_ps(&m_ReversedSignal[i]);
			kernel = _mm_load_ps(&m_Taps[i]);
			product = _mm_mul_ps(signal, kernel);
			accumulator = _mm_add_ps(product, accumulator);
		}

		output = accumulator.m128_f32[0] +
			accumulator.m128_f32[1] +
			accumulator.m128_f32[2] +
			accumulator.m128_f32[3];

		// finish the tail:
		for (int j = MultOf4Size; j < size; j++) {
			output += m_ReversedSignal[j] * m_Taps[j];
		}

#else

		// unroll 2 (fastest):	

		int index = m_CurrentIndex;
		int i;

		for (i = 0; i < (size >> 1) << 1; i += 2) {

			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_Taps[i];

			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_Taps[i + 1];
		}

		// Tail:
		if (size & 1) { // Do one more if odd number:
			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_Taps[i];
		}

		// unroll 4:
		//int index = m_CurrentIndex;
		//int i;
		//
		//int MultOf4Size = (size >> 2) << 2;
		//for (i = 0; i < MultOf4Size; i += 4) {
		//	index = (index == 0) ? size - 1 : index - 1;
		//	output += m_Signal[index] * m_Taps[i];
		//	index = (index == 0) ? size - 1 : index - 1;
		//	output += m_Signal[index] * m_Taps[i + 1];
		//	index = (index == 0) ? size - 1 : index - 1;
		//	output += m_Signal[index] * m_Taps[i + 2];
		//	index = (index == 0) ? size - 1 : index - 1;
		//	output += m_Signal[index] * m_Taps[i + 3];
		//}
		//// Tail:
		//for (int j = MultOf4Size; j < size; j++) {
		//	index = (index == 0) ? size - 1 : index - 1;
		//	output += m_Signal[index] * m_Taps[i];
		//}

#endif
		return output;
	}

private:
	alignas(16) FloatType m_Taps[size];
	alignas(16) FloatType m_Signal[size];
	int JobLength;
	int JobTailLength;
	int m_CurrentIndex;

#ifdef USE_SIMD
	alignas(16) FloatType m_ReversedSignal[size];
#endif
	
};


#endif