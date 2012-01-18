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

#ifndef LATL_SE2_HPP
#define LATL_SE2_HPP

#include <latl/so2.hpp>

namespace latl
{    
    template <class Scalar=double>
    class SE2 {
    private:
        Vector<2,Scalar> t;
        SO2<Scalar> R;

        static void compute_coef(Scalar theta,
                                 Scalar& A, Scalar& B)
        {
            Scalar theta_sq = theta * theta;
            if (theta_sq < Constants<Scalar>::sqrt_epsilon()) {
                A = 1 - theta_sq*Scalar(1.0/6) * (1 + theta_sq*Scalar(0.05));
                B = theta*Scalar(0.5)*(1 - theta_sq*Scalar(1.0/12)*(1 - theta_sq*Scalar(1.0/30)));
            } else {
                Scalar inv_theta = 1/theta;
                A = latl::sin(theta) * inv_theta;
                B = (1 - latl::cos(theta)) * inv_theta;
            }
        }
        
    public:
        typedef Scalar scalar_type;
        enum {DoF = 3};
        
        SE2() : t(Zeros()) {}

        SE2(const Uninitialized& u) {}
        
        template <class V, class S>
        SE2(const AbstractVector<V>& t_, const SO2<S>& R_)
            : t(t_), R(R_) {}

        template <class V>
        SE2(const AbstractVector<V>& t_)
            : t(t_){}

        SE2(const SO2<Scalar>& R_)
            : t(Zeros()), R(R_){}

        SE2(const SE2& other) : t(other.t), R(other.R) {}

        template <class T>
        explicit SE2(const SE2<T>& other) : t(other.t), R(other.R) {}
        
        void rectify() {
            R.rectify();
        }

        const Vector<2,Scalar> translation() const { return t; }
        Vector<2,Scalar>& translation() { return t; }

        const SO2<Scalar>& rotation() const { return R; }
        SO2<Scalar>& rotation() { return R; }
        
        template <class V>
        Vector<2,Scalar> operator*(const FixedVector<2,V>& x) const {
            return R * x + t;
        }

        template <class V>
        Vector<3,Scalar> operator*(const FixedVector<3,V>& x) const {
            Vector<3,Scalar> y;
            slice<0,2>(y) = R * slice<0,2>(x) + t * x[2];
            y[2] = x[2];
            return y;
        }
                
        SE2 operator*(const SE2& other) const {
            return SE2(t + R*other.t, R * other.R);
        }

        SE2 inverse() const {
            return SE2(-(R.inverse()*t), R.inverse());
        }
        
        Vector<3,Scalar> ln() const {
            Scalar theta = R.ln();
            Scalar A,B;
            compute_coef(theta, A, B);

            Scalar f = 1 / (A*A + B*B);
            Vector<3,Scalar> uw;
            uw[0] = f * ( A * t[0] + B * t[1]);
            uw[1] = f * (-B * t[0] + A * t[1]);
            uw[2] = theta;
            return uw;
        }

        template <class V>
        static SE2 exp(const AbstractVector<V>& uw) {
            assert_size_is<3>(uw);
            Scalar A,B;
            compute_coef(uw[2], A, B);

            SE2<Scalar> exp_uw((Uninitialized()));
            SO2<Scalar>::compute_exp_matrix(uw[2], exp_uw.R.R);
            
            exp_uw.translation()[0] = A*uw[0] - B*uw[1];
            exp_uw.translation()[1] = B*uw[0] + A*uw[1];
            
            return exp_uw;
        }

        template <class Mat>
        void adjoint(AbstractMatrix<Mat>& adj) const
        {
            assert_shape_is<3,3>(adj);
            slice<0,0,2,2>(adj) = R.matrix();
            adj(0,2) =  t[1];
            adj(1,2) = -t[0];
            adj(2,0) = adj(2,1) = 0;
            adj(2,2) = 1;
        }
        
        template <class V>
        Vector<3,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<3>(v);
            Vector<3,Scalar> adj_v;
            slice<0,2>(adj_v) = R * slice<0,2>(v);
            adj_v[0] += t[1] * v[2];
            adj_v[1] -= t[0] * v[2];
            adj_v[2] = v[2];
            return adj_v;
        }

        template <class V>
        Vector<3,Scalar> adjoint_T_times(const AbstractVector<V>& v) const {
            assert_size_is<3>(v);
            Vector<3,Scalar> adjT_v;
            slice<0,2>(adjT_v) = R.matrix().T() * slice<0,2>(v);
            adjT_v[2] = t[1] * v[0] - t[0] * v[1] + v[2];
            return adjT_v;
        }
    };
}

#endif
