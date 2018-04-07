ReSampler-QuadMath.exe
======================

Experimental Quadruple-precision version of ReSampler

* uses Quadruple-Precision (128-bit floating-point) for calculations
* very slow, because CPU doesn't have 128-bit floating point, so it's implemented in software
* uses the GNU quad-precision math library - requires libquadmath-0.dll
* Currently no method of reading quad-precision input files. (The best you can do is 64f)

Additional Dependencies
=======================
libquadmath-0.dll