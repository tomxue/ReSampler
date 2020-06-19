#ifndef IQDEMODULATOR_H
#define IQDEMODULATOR_H

template<typename FloatType>
class IQDemodulator
{

public:
    FloatType demodulateFM(FloatType i, FloatType q)
    {
        i2 = i1;
        i1 = i0;
        i0 = i;
        q2 = q1;
        q1 = q0;
        q0 = q;

        return ((q0 - q2) * i1) - ((i0 - i2) * q1);
    }

private:
    FloatType i0{0.0};
    FloatType i1{0.0};
    FloatType i2{0.0};
    FloatType q0{0.0};
    FloatType q1{0.0};
    FloatType q2{0.0};

};

#endif // IQDEMODULATOR_H
