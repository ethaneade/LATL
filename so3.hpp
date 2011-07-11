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

#ifndef LATL_SO3_HPP
#define LATL_SO3_HPP

#include <latl/latl.hpp>

namespace latl
{
    template <class Scalar>
    class SE3;
    
    template <class Scalar>
    class SO3 {
    private:
        Matrix<3,3,Scalar> R;
        explicit SO3(const Matrix<3,3,Scalar>& R_) : R(R_) {}
        friend class SE3<Scalar>;
    public:
        SO3() : R(Identity()) {}

        SO3(const SO3& other) : R(other.R) {}

        template <class T>
        explicit SO3(const SO3<T>& other) : R(other.R) {}
        
        template <class Mat>
        static SO3 from_matrix(const AbstractMatrix<Mat>& m) {
            Matrix<3,3,Scalar> R = m;
            orthonormalize(R);
            return SO3(R);
        }
        
        void rectify() {
            orthonormalize(R);
        }
        
        template <class V>
        Vector<3,Scalar> operator*(const AbstractVector<V>& x) const {
            return R * x;
        }

        SO3 operator*(const SO3& other) const {
            return SO3(R * other.R);
        }

        SO3 inverse() const {
            return SO3(R.T());
        }

        const Matrix<3,3,Scalar>& matrix() const { return R; }

        Vector<3,Scalar> ln() const {
            Vector<3,Scalar> w;
            w[0] = R(2,1) - R(1,2);
            w[1] = R(0,2) - R(2,0);
            w[2] = R(1,0) - R(0,1);
            
            Scalar cos_theta = latl::max(Scalar(-1),
                                         latl::min(Scalar(0.5)*(trace(R) - 1),
                                                   Scalar(1)));
            Scalar theta = latl::acos(cos_theta);
            
            if (cos_theta > Scalar(-0.99)) {
                return w / (2*sinc(theta));
            }

            Scalar inv_B = (theta * theta) / (1 - cos_theta);
            
            w[0] = latl::sqrt(inv_B * (R(0,0) - cos_theta));
            w[1] = latl::sqrt(inv_B * (R(1,1) - cos_theta));
            w[2] = latl::sqrt(inv_B * (R(2,2) - cos_theta));

            if (R(2,1) < R(1,2))
                w[0] = -w[0];
            if (R(0,2) < R(2,0))
                w[1] = -w[1];
            if (R(1,0) < R(0,1))
                w[2] = -w[2];

            return w;
        }

        template <class V>
        static SO3 exp(const AbstractVector<V>& w) {
            assert_size_is<3>(w);
            Scalar theta_sq = norm_sq(w);
            Scalar A,B;
            if (theta_sq < Constants<Scalar>::sqrt_epsilon()) {
                A = 1 - theta_sq*Scalar(1.0/6) * (1 + theta_sq*Scalar(0.05));
                B = Scalar(0.5)*(1 - theta_sq*Scalar(1.0/12)*(1 - theta_sq*Scalar(1.0/30)));
            } else {
                Scalar theta = latl::sqrt(theta_sq);
                A = latl::sin(theta) / theta;
                B = (1 - latl::cos(theta)) / theta_sq;
            }
            
            Matrix<3,3,Scalar> wx = skew_matrix(w);
            Matrix<3,3,Scalar> R = Identity();
            R += A*wx;
            R += B*(wx*wx);
            return SO3(R);
        }

        template <class V>
        Vector<3,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<3>(v);
            return R * v;
        }

        template <class V>
        Vector<3,Scalar> adjoint_T_times(const AbstractVector<V>& v) const {
            assert_size_is<3>(v);
            return R.T() * v;
        }        
    };
}

#endif
