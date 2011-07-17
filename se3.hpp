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

#ifndef LATL_SE3_HPP
#define LATL_SE3_HPP

#include <latl/so3.hpp>
#include <latl/lu.hpp>

namespace latl
{    
    template <class Scalar>
    class SE3 {
    private:
        Vector<3,Scalar> t;
        SO3<Scalar> R;

        static void compute_coef(Scalar theta_sq,
                                 Scalar& A, Scalar& B, Scalar& C)
        {
            if (theta_sq < Constants<Scalar>::sqrt_epsilon()) {
                A = 1 - theta_sq*Scalar(1.0/6) * (1 + theta_sq*Scalar(0.05));
                B = Scalar(0.5)*(1 - theta_sq*Scalar(1.0/12)*(1 - theta_sq*Scalar(1.0/30)));
                C = Scalar(1.0/6)*(1 - theta_sq*Scalar(0.05)*(1 - theta_sq*Scalar(1.0/42)));
            } else {
                Scalar theta = latl::sqrt(theta_sq);
                Scalar inv_theta_sq = 1/theta_sq;
                A = latl::sin(theta) / theta;
                B = (1 - latl::cos(theta)) * inv_theta_sq;
                C = (1 - A) * inv_theta_sq;
            }            
        }
        
    public:
        SE3() : t(Zeros()) {}

        SE3(const Uninitialized& u) {}
        
        template <class V, class S>
        SE3(const AbstractVector<V>& t_, const SO3<S>& R_)
            : t(t_), R(R_) {}

        template <class V>
        SE3(const AbstractVector<V>& t_)
            : t(t_){}

        SE3(const SO3<Scalar>& R_)
            : t(Zeros()), R(R_){}

        SE3(const SE3& other) : t(other.t), R(other.R) {}

        template <class T>
        explicit SE3(const SE3<T>& other) : t(other.t), R(other.R) {}
        
        void rectify() {
            R.rectify();
        }

        const Vector<3,Scalar> translation() const { return t; }
        Vector<3,Scalar>& translation() { return t; }

        const SO3<Scalar>& rotation() const { return R; }
        SO3<Scalar>& rotation() { return R; }
        
        template <class V>
        Vector<3,Scalar> operator*(const FixedVector<3,V>& x) const {
            return R * x + t;
        }

        template <class V>
        Vector<4,Scalar> operator*(const FixedVector<4,V>& x) const {
            Vector<4,Scalar> y;
            slice<0,3>(y) = R * slice<0,3>(x) + t * x[3];
            y[3] = x[3];
            return y;
        }
                
        SE3 operator*(const SE3& other) const {
            return SE3(t + R*other.t, R * other.R);
        }

        SE3 inverse() const {
            return SE3(-(R.inverse()*t), R.inverse());
        }
        
        Vector<6,Scalar> ln() const {
            Vector<3,Scalar> w = R.ln();
            Scalar theta_sq = norm_sq(w);
            Scalar A,B,C;
            compute_coef(theta_sq, A, B, C);
            
            Matrix<3,3,Scalar> Q;
            
            SO3<Scalar>::compute_exp_matrix(w, B, C, Q);
            LU<3,Scalar> luQ(Q);
            Vector<3,Scalar> u = luQ.inverse_times(t);
            return (u, w);
        }

        template <class V>
        static SE3 exp(const AbstractVector<V>& uw) {
            assert_size_is<6>(uw);
            Vector<3,Scalar> u = slice<0,3>(uw);
            Vector<3,Scalar> w = slice<3,3>(uw);
            Scalar theta_sq = norm_sq(w);
            Scalar A,B,C;
            compute_coef(theta_sq, A, B, C);

            SE3<Scalar> exp_uw((Uninitialized()));
            SO3<Scalar>::compute_exp_matrix(w, A, B, exp_uw.R.R);
            
            Vector<3,Scalar> wxu = w ^ u;
            exp_uw.translation() = u +  B * wxu + C * (w ^ wxu);
            return exp_uw;
        }

        template <class V>
        Vector<6,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<6>(v);
            Vector<3,Scalar> Rw = R * slice<3,3>(v);
            return (R * slice<0,3>(v) + (t ^ Rw), Rw);
        }

        template <class V>
        Vector<6,Scalar> adjoint_T_times(const AbstractVector<V>& v) const {
            assert_size_is<6>(v);
            Vector<3,Scalar> RTu = R.T() * slice<0,3>(v);            
            return (R.T() * slice<0,3>(v),
                    R.T() * (slice<3,3>(v) - t ^ slice<0,3>(v)));         
        }
    };

    template <class Scalar, class Mat>
    void transform_covariance(const SE3<Scalar>& trans,
                              AbstractMatrix<Mat>& m)
    {
        assert_square(m);
        assert_size_is<6>(m[0]);

        for (int i=0; i<6; ++i)
            m.T()[i] = trans.adjoint_times(m.T()[i]);
        for (int i=0; i<6; ++i)
            m[i] = trans.adjoint_times(m[i]);
    }

    template <class Scalar, class Mat>
    void transform_information(const SE3<Scalar>& trans,
                               AbstractMatrix<Mat>& m)
    {
        assert_square(m);
        assert_size_is<6>(m[0]);

        for (int i=0; i<6; ++i)
            m.T()[i] = trans.adjoint_T_times(m.T()[i]);
        for (int i=0; i<6; ++i)
            m[i] = trans.adjoint_T_times(m[i]);
    }
    
}

#endif
