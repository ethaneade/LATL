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

#ifndef LATL_SIM2_HPP
#define LATL_SIM2_HPP

#include <latl/se2.hpp>

namespace latl
{    
    template <class Scalar=double>
    class Sim2 {
    private:
        SE2<Scalar> rigid_;
        Scalar scale_;

        static void compute_coef(Scalar theta,
                                 Scalar& A, Scalar& B, Scalar& C)
        {
            Scalar theta_sq = theta * theta;
            if (theta_sq < Constants<Scalar>::sqrt_epsilon()) {
                A = 1 - theta_sq*Scalar(1.0/6) * (1 + theta_sq*Scalar(0.05));
                B = Scalar(0.5)*(1 - theta_sq*Scalar(1.0/12)*(1 - theta_sq*Scalar(1.0/30)));
                C = Scalar(1.0/6)*(1 - theta_sq*Scalar(0.05)*(1 - theta_sq*Scalar(1.0/42)));
            } else {
                Scalar inv_theta = 1/theta;
                Scalar inv_theta_sq = inv_theta * inv_theta;
                A = latl::sin(theta) * inv_theta;
                B = (1 - latl::cos(theta)) * inv_theta_sq;
                C = (1 - A) * inv_theta_sq;
            }
        }

        static void compute_V(Scalar theta, Scalar lambda, Scalar exp_lambda,
                              Matrix<2,2,Scalar>& V)
        {
            Scalar lambda_sq = lambda*lambda;            
            Scalar alpha_den = lambda_sq + theta*theta;
            Scalar alpha;
            if (alpha_den == 0)
                alpha = 1;
            else
                alpha = lambda_sq / alpha_den;
            
            Scalar A,B,C;
            compute_coef(theta, A, B, C);

            Scalar exp_neg_lambda = 1/exp_lambda;
                        
            Scalar L;
            Scalar beta;
            if (lambda_sq < Constants<Scalar>::sqrt_epsilon()) {
                L = 1 - lambda*Scalar(0.5)*(1 - lambda*Scalar(1.0/3)*(1 - lambda*Scalar(0.25)*(1 - lambda*Scalar(0.2))));
                beta = 0.5*(1 - lambda*Scalar(1.0/3)*(1 - lambda*Scalar(0.25)*(1 - lambda*Scalar(0.2)*(1 - lambda*Scalar(1.0/6)))));
            } else {
                L = (1 - exp_neg_lambda) / lambda;
                beta = (exp_neg_lambda - 1 + lambda) / lambda_sq;
            }

            Scalar z = A - lambda*B;
            Scalar gamma = B - lambda*C;
            Scalar a = alpha*(L - z) + z;
            Scalar c = theta * (alpha*(beta - gamma) + gamma);

            V(0,0) = V(1,1) = a;
            V(0,1) = -c;
            V(1,0) = c;
        }
        
    public:
        typedef Scalar scalar_type;
        enum {DoF = 4};
        
        Sim2(Scalar scale=1) : scale_(scale) {}

        Sim2(const Uninitialized& u) : rigid_(u) {}
        
        Sim2(const SE2<Scalar>& rigid, Scalar s=1)
            : rigid_(rigid), scale_(s) {}
        
        template <class T>
        explicit Sim2(const Sim2<T>& other) : rigid_(other.rigid_)
        {
            scale_ = (Scalar)other.scale_;
        }
        
        void rectify() {
            rigid_.rectify();
        }

        SE2<Scalar>& rigid() { return rigid_; }
        const SE2<Scalar>& rigid() const { return rigid_; }
        
        const Vector<2,Scalar> translation() const { return rigid_.translation(); }
        Vector<2,Scalar>& translation() { return rigid_.translation(); }

        const SO2<Scalar>& rotation() const { return rigid_.rotation(); }
        SO2<Scalar>& rotation() { return rigid_.rotation(); }

        Scalar scale() const { return scale_; }
        Scalar& scale() { return scale_; }
        
        template <class V>
        Vector<2,Scalar> operator*(const FixedVector<2,V>& x) const {
            return scale_ * (rigid_ * x);
        }

        template <class V>
        Vector<3,Scalar> operator*(const FixedVector<3,V>& x) const {
            Vector<3,Scalar> y = rigid_ * x;
            slice<0,2>(y) *= scale_;
            return y;
        }
                
        Sim2 operator*(const Sim2& other) const {
            Sim2 p((Uninitialized()));
            p.rigid_.rotation() = rotation() * other.rotation();
            p.rigid_.translation() = rotation() * other.translation();
            p.rigid_.translation() += translation() / other.scale();
            p.scale() = scale() * other.scale();
            return p;
        }

        Sim2 inverse() const {            
            return Sim2(SE2<Scalar>(-scale_*(rotation().matrix().T()*translation()),
                                    rotation().inverse()),
                        1/scale_);
        }
        
        Vector<4,Scalar> ln() const {
            Scalar theta = rotation().ln();
            Scalar lambda = latl::ln(scale_);

            Matrix<2,2,Scalar> V;
            compute_V(theta, lambda, scale_, V);

            Vector<4,Scalar> uwl;
            slice<0,2>(uwl) = latl::inverse(V) * translation();
            uwl[2] = theta;
            uwl[3] = lambda;
            
            return uwl;
        }

        template <class V>
        static Sim2 exp(const AbstractVector<V>& uwl) {
            assert_size_is<4>(uwl);
            Scalar theta = uwl[2];
            Scalar lambda = uwl[3];

            Matrix<2,2,Scalar> V;
            Scalar exp_lambda = latl::exp(lambda);
            compute_V(theta, lambda, exp_lambda, V);
            
            Sim2<Scalar> exp_uw((Uninitialized()));
            exp_uw.rotation() = SO2<Scalar>::exp(theta);
            exp_uw.scale() = exp_lambda;            
            exp_uw.translation() = V * slice<0,2>(uwl);
            
            return exp_uw;
        }

        template <class Mat>
        void adjoint(AbstractMatrix<Mat>& adj) const
        {
            assert_shape_is<4,4>(adj);
            slice<0,0,2,2>(adj) = scale_ * rotation().matrix();
            const Vector<2,Scalar>& t = translation();
            adj(0,2) = scale_ * t[1];
            adj(0,3) = scale_ *-t[0];
            adj(1,2) = scale_ *-t[0];
            adj(1,3) = scale_ *-t[1];
            adj(2,2) = adj(3,3) = 1;
            adj(2,3) = adj(3,2) = 0;
            zero(slice<2,0,2,2>(adj).instance());
        }
        
        template <class V>
        Vector<4,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<4>(v);
            Vector<4,Scalar> adj_v;
            const Vector<2,Scalar>& t = translation();

            slice<0,2>(adj_v) = scale_ * (rotation() * slice<0,2>(v));
            adj_v[0] += scale_ * ( t[1] * v[2] - t[0] * v[3]);
            adj_v[1] += scale_ * (-t[0] * v[2] - t[1] * v[3]);
            adj_v[2] = v[2];
            adj_v[3] = v[3];
            return adj_v;
        }

        template <class V>
        Vector<4,Scalar> adjoint_T_times(const AbstractVector<V>& v) const {
            assert_size_is<4>(v);
            Vector<4,Scalar> adjT_v;
            const Vector<2,Scalar>& t = translation();
            slice<0,2>(adjT_v) = scale_*(rotation().matrix().T() * slice<0,2>(v));
            adjT_v[2] = v[2] + scale_*(t[1]*v[0] - t[0]*v[1]);
            adjT_v[3] = v[3] - scale_*(t[0]*v[0] + t[1]*v[1]);
            
            return adjT_v;
        }
    };
}

#endif
