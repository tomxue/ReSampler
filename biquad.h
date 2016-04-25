#ifndef biquad_H
#define biquad_H

// Simple iir biquad filter implementation

template<typename FloatType> class Biquad {
public:
	
	// Naming conventions for coefficients are as follows:
	// numerator/input coeffs:			a0,a1,a2
	// denominator/feedback coeffs:		b1,b2
	Biquad() {};
	Biquad(FloatType a0, FloatType a1, FloatType a2, FloatType b1, FloatType b2)
		: a0(a0), a1(a1), a2(a2), b1(b1), b2(b2),z1(0.0),z2(0.0) {};

	FloatType filter(FloatType inSample) {

		// Biquad calculation using transposed Direct Form 2:
		
		FloatType outSample = inSample * a0 + z1;
		z1 = inSample * a1 + z2 - b1 * outSample;
		z2 = inSample * a2 - b2 * outSample;
		return outSample;
	}

	void setCoeffs(FloatType a0, FloatType a1, FloatType a2, FloatType b1, FloatType b2) {
		Biquad::a0 = a0;
		Biquad::a1 = a1;
		Biquad::a2 = a2;
		Biquad::b1 = b1;
		Biquad::b2 = b2;
	}

    
protected:
    FloatType a0, a1, a2, b1, b2;
    double z1, z2;
};

#endif biquad_H
