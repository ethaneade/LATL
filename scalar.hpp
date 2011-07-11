// Copyright 2011 Ethan Eade. All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:

//    1. Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.

//    2. Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY ETHAN EADE ``AS IS'' AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ETHAN EADE OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

// The views and conclusions contained in the software and
// documentation are those of the authors and should not be
// interpreted as representing official policies, either expressed or
// implied, of Ethan Eade.

#ifndef LATL_SCALAR_HPP
#define LATL_SCALAR_HPP

#include <latl/debug.hpp>
#include <float.h>
#include <cmath>

namespace latl
{

    template <class T>
    struct ScalarType;

    template <> struct
    ScalarType<int> {
        typedef int type;
        typedef int nonrecip_type;
    };

    template <> struct
    ScalarType<float> {
        typedef float type;
        typedef float recip_type;
    };

    template <> struct
    ScalarType<double> {
        typedef double type;
        typedef double recip_type;
    };
    
    template <class A, class B>
    struct Wider;

    template <class A>
    struct Wider<A,A> { typedef A type; };

    template <> struct Wider<int,float> { typedef float type; };
    template <> struct Wider<int,double> { typedef double type; };
        
    template <> struct Wider<float,int> { typedef float type; };
    template <> struct Wider<float,double> { typedef double type; };

    template <> struct Wider<double,int> { typedef double type; };
    template <> struct Wider<double,float> { typedef double type; };

    template <class T> typename ScalarType<T>::type
    min(T a, T b) { return a < b ? a : b; }
    
    template <class T> typename ScalarType<T>::type
    max(T a, T b) { return a > b ? a : b; }
    
    template <class Scalar>
    struct Constants;

    template <class Scalar>
    struct ScalarFuncs;

    template <>
    struct Constants<float> {
        static float epsilon() { return 1.192093e-07f; }        
        static float sqrt_epsilon() { return 3.452670e-04f; }
        static int digits() { return 7; }
    };

    template <>
    struct Constants<double> {
        static double epsilon() { return 2.220446e-16; }
        static double sqrt_epsilon() { return 1.490116e-08; }
        static int digits() { return 16; }
    };
    
    template <>
    struct ScalarFuncs<float> {
        static float abs(float s) { return fabsf(s); }
        static float sqrt(float s) { return sqrtf(s); }
        static float sin(float s) { return sinf(s); }
        static float cos(float s) { return cosf(s); }
        static float acos(float s) { return acosf(s); }
        static float tan(float s) { return tanf(s); }
        static float atan(float s) { return atanf(s); }
        static float atan2(float y, float x) { return atan2f(y,x); }
    };
        
    template <>
    struct ScalarFuncs<double> {
        static double abs(double s) { return ::fabs(s); }
        static double sqrt(double s) { return ::sqrt(s); }
        static double sin(double s) { return ::sin(s); }
        static double cos(double s) { return ::cos(s); }
        static double acos(double s) { return ::acos(s); }
        static double tan(double s) { return ::tan(s); }
        static double atan(double s) { return ::atan(s); }
        static double atan2(double y, double x) { return ::atan2(y,x); }
    };

    template <class T> typename ScalarType<T>::type
    abs(T x) { return ScalarFuncs<T>::abs(x); }

    template <class T> typename ScalarType<T>::type
    sqrt(T x) { return ScalarFuncs<T>::sqrt(x); }

    template <class T> typename ScalarType<T>::type
    sin(T x) { return ScalarFuncs<T>::sin(x); }

    template <class T> typename ScalarType<T>::type
    cos(T x) { return ScalarFuncs<T>::cos(x); }

    template <class T> typename ScalarType<T>::type
    acos(T x) { return ScalarFuncs<T>::acos(x); }
    
    template <class T> typename ScalarType<T>::type
    tan(T x) { return ScalarFuncs<T>::tan(x); }

    template <class T> typename ScalarType<T>::type
    atan(T x) { return ScalarFuncs<T>::atan(x); }
    
    template <class T> typename ScalarType<T>::type
    atan2(T y, T x) { return ScalarFuncs<T>::atan2(y, x); }

    template <class T> typename ScalarType<T>::type
    sq(T x) { return x*x; }

    template <class T> typename ScalarType<T>::type
    sinc(T x) {
        T xx = sq(x);
        T x4 = sq(xx);
        if (x4 < Constants<T>::epsilon())
            return 1 - xx*T(1.0/6) * (1 + xx*T(0.05));
        else
            return latl::sin(x) / x;
    }
}

#endif
