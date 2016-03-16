#ifndef FIRFILTER_H_
#define FIRFFILTER_H_

// #define VECTORIZE_FIR_GET 1 13/03/2016; Still NQR !

template <typename FloatType, unsigned int size>
class FirFilter {
public:

	FloatType* m_pTaps;

	FirFilter(FloatType* taps) :
		m_CurrentIndex(0) {
		m_pTaps = taps;
		for (unsigned int i = 0; i < size; ++i)
			m_Signal[i] = 0.0;
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

#ifdef VECTORIZE_FIR_GET


		int IndexA = (m_CurrentIndex);
		int IndexB = (m_CurrentIndex - 1);
		int IndexC = (m_CurrentIndex - 2);
		int IndexD = (m_CurrentIndex - 3);


		FloatType accA(0), accB(0), accC(0), accD(0); // accumlators

		for (int i = 0; i < size >> 2; ++i) {

			// 32766 size-1
			// 32765 size-2													/
			// 32764 size-3													// 0->size-5, -1->size-6, -2->size-7 -3->size-8 
			// 32763 size-4													// 4->size-1, 3->size-2, 2->size-3, 1->Size-4 

			IndexA = ((IndexA - 4) <= 0) ? size - 5 + IndexA : IndexA - 4;
			IndexB = ((IndexB - 4) <= 0) ? size - 5 + IndexB : IndexB - 4;
			IndexC = ((IndexC - 4) <= 0) ? size - 5 + IndexC : IndexC - 4;
			IndexD = ((IndexD - 4) <= 0) ? size - 5 + IndexD : IndexD - 4;

			accA += m_Signal[IndexA] * m_pTaps[i << 2 /* i*4 */];
			accB += m_Signal[IndexB] * m_pTaps[1 + (i << 2) /* 1+i*4 */];
			accC += m_Signal[IndexC] * m_pTaps[2 + (i << 2) /* 2+i*4 */];
			accD += m_Signal[IndexD] * m_pTaps[3 + (i << 2) /* 3+i*4 */];

		}
		// finish tail
		for (int j = (size >> 2) << 2; j < size; j++)
		{
			IndexD = (IndexD == 0) ? size - 1 : IndexD - 1;
			accD += m_Signal[IndexD] * m_pTaps[j];
		}

		output = accA + accB + accC + accD;

#else
		int index = m_CurrentIndex;
		int i;
		// unroll 2:
		for (i = 0; i < (size >> 1) << 1; i += 2) {

			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_pTaps[i];

			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_pTaps[i + 1];
		}

		// Tail:
		if (size & 1) { // Do one more if odd number:
			index = (index == 0) ? size - 1 : index - 1;
			output += m_Signal[index] * m_pTaps[i];
		}
#endif
		return output;
	}

private:
	FloatType m_Signal[size];
	int m_CurrentIndex;
};

#endif