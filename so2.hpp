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

#ifndef LATL_SO2_HPP
#define LATL_SO2_HPP

#include <latl/latl.hpp>
#include <latl/uninitialized.hpp>

namespace latl
{
    template <class Scalar>
    class SE2;
    
    template <class Scalar>
    class SO2 {
    private:
        Matrix<2,2,Scalar> R;
        explicit SO2(const Matrix<2,2,Scalar>& R_) : R(R_) {}
        friend class SE2<Scalar>;

        static void compute_exp_matrix(const Scalar theta,
                                       Matrix<2,2,Scalar>& m)
        {
            Scalar st = latl::sin(theta);
            Scalar ct = latl::cos(theta);

            m(0,0) = m(1,1) = ct;
            m(0,1) = -st;
            m(1,0) = st;
        }
        
    public: 
        SO2() : R(Identity()) {}

        SO2(const Uninitialized& u) {}

        SO2(const SO2& other) : R(other.R) {}

        template <class T>
        explicit SO2(const SO2<T>& other) : R(other.R) {}
        
        template <class Mat>
        static SO2 from_matrix(const AbstractMatrix<Mat>& m) {
            Matrix<2,2,Scalar> R = m;
            orthonormalize(R);
            return SO2(R);
        }

        template <class Mat>
        static SO2 bless(const AbstractMatrix<Mat>& m) {
            return SO2(m);
        }
        
        void rectify() {
            orthonormalize(R);
        }
        
        template <class V>
        Vector<2,Scalar> operator*(const AbstractVector<V>& x) const {
            return R * x;
        }

        SO2 operator*(const SO2& other) const {
            return SO2(R * other.R);
        }

        SO2 inverse() const {
            return SO2(R.T());
        }

        const Matrix<2,2,Scalar>& matrix() const { return R; }

        Scalar ln() const {
            return latl::atan2(R(1,0), R(0,0));
        }

        static SO2 exp(Scalar theta) {
            SO2 expw((Uninitialized()));
            compute_exp_matrix(theta, expw.R);
            return expw;
        }
    };
}

#endif
