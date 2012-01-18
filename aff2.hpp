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

#ifndef LATL_AFF2_HPP
#define LATL_AFF2_HPP

#include <latl/sim2.hpp>

namespace latl
{    
    template <class Scalar=double>
    class Aff2 {
    private:
        Matrix<2,2,Scalar> A;
        Vector<2,Scalar> t;
                
    public:
        typedef Scalar scalar_type;
        enum {DoF = 6};
        
        Aff2() : A(Identity()), t(Scalar(0)) {}

        Aff2(const Uninitialized& u) {}

        template <class Mat>
        Aff2(const FixedMatrix<2,3,Mat>& At)
            : A(slice<0,0,2,2>(At)), t(At.T()[2]) {}
        
        template <class Mat, class V>
        Aff2(const FixedMatrix<2,2,Mat>& A_, const FixedVector<2,V>& t_)
            : A(A_), t(t_) {}

        template <class V>
        Aff2(const FixedVector<2,V>& t_)
            : A(Identity()), t(t_) {}

        template <class Mat>
        Aff2(const FixedMatrix<2,2,Mat>& A_)
            : A(A_), t(Zeros()) {}
        
        Aff2(const SE2<Scalar>& se2)
            : A(se2.rotation().matrix()), t(se2.translation()) {}

        Aff2(const Sim2<Scalar>& sim2)
            : A(sim2.scale() * sim2.rotation().matrix()),
              t(sim2.scale() * sim2.translation()) {}
                
        template <class T>
        explicit Aff2(const Aff2<T>& other) : A(other.A), t(other.t) {}
        
        const Vector<2,Scalar> translation() const { return t; }
        Vector<2,Scalar>& translation() { return t; }

        const Matrix<2,2,Scalar>& linear() const { return A; }
        Matrix<2,2,Scalar>& linear() { return A; }

        template <class V>
        Vector<2,Scalar> operator*(const FixedVector<2,V>& x) const {
            return A*x + t;
        }

        template <class V>
        Vector<3,Scalar> operator*(const FixedVector<3,V>& x) const {
            Vector<3,Scalar> y;
            slice<0,2>(y) = A*slice<0,2>(x) + x[2]*t;
            y[2] = x[2];
            return y;
        }
                
        Aff2 operator*(const Aff2& other) const {
            return Aff2(A*other.A, A*other.t + t);
        }
        
        Aff2 inverse() const {
            Matrix<2,2,Scalar> invA = latl::inverse(A);
            return Aff2(invA, -(invA*t));
        }

        template <class V>
        bool ln(FixedVector<6,V>& v) const {
            Matrix<3,3,Scalar> M;
            slice<0,0,2,2>(M) = A;
            slice<0,2>(M.T()[2].instance()) = t;
            M(2,0) = M(2,1) = 0;
            M(2,2) = 1;
            if (!log(M, M))
                return false;
            slice<0,2>(v) = slice<0,2>(M.T()[2]);
            v[2] = Scalar(0.5)*(M(1,0) - M(0,1));
            v[3] = Scalar(0.5)*(M(0,0) + M(1,1));
            v[4] = Scalar(0.5)*(M(0,0) - M(1,1));
            v[5] = Scalar(0.5)*(M(1,0) + M(0,1));
            return true;
        }

        Vector<6,Scalar> ln() const {
            Vector<6,Scalar> v;
            bool success = ln(v);
            assert(success);
            return v;
        }
        
        template <class V>
        static Aff2 exp(const AbstractVector<V>& v) {
            assert_size_is<6>(v);
            Matrix<3,3,Scalar> M;
            M(0,2) = v[0];
            M(1,2) = v[1];
            M(0,0) = v[3] + v[4];
            M(0,1) = v[5] - v[2];
            M(1,0) = v[5] + v[2];
            M(1,1) = v[3] - v[4];
            M(2,0) = M(2,1) = M(2,2) = 0;
            latl::exp(M, M);            
            
            return Aff2<Scalar>(slice<0,0,2,3>(M));
        }

        template <class Mat>
        void adjoint(AbstractMatrix<Mat>& adj) const
        {
            assert_shape_is<6,6>(adj);
            Scalar s = 1/determinant(A);
            Scalar a = A(0,0), b = A(0,1), c = A(1,0), d = A(1,1);
            Scalar aa=a*a, bb=b*b, cc=c*c, dd=d*d;
            Scalar ab=a*b, ac=a*c, ad=a*d;
            Scalar bc=b*c, bd=b*d, cd=c*d;

            Scalar half_s = s * Scalar(0.5);
            adj(2,2) = half_s*( aa+bb+cc+dd);
            adj(2,5) = half_s*(-aa+bb-cc+dd);
            adj(5,2) = half_s*(-aa-bb+cc+dd);
            adj(5,5) = half_s*( aa-bb-cc+dd);

            adj(2,3) = 0;
            adj(3,2) = adj(3,4) = adj(3,5) = 0;
            adj(4,3) = adj(5,3) = 0;

            Scalar ab_cd = ab+cd;
            Scalar ac_bd = ac+bd;
            Scalar ad_bc = ad+bc;
            
            adj(2,4) = s*(ab_cd);
            adj(3,3) = 1;
            adj(4,2) = s*(ac_bd);
            adj(4,4) = s*(ad_bc);
            adj(4,5) = s*(-ac+bd);
            adj(5,4) = s*(-ab+cd);
            
            slice<0,0,2,2>(adj) = A;

            adj(0,2) = s*((aa+bb)*t[1] - ac_bd*t[0]);
            adj(0,3) = -t[0];
            adj(0,4) = s*(2*ab*t[1] - ad_bc*t[0]);
            adj(0,5) = s*((ac-bd)*t[0] - (aa-bb)*t[1]);

            adj(1,2) = s*((ac_bd)*t[1] - (cc+dd)*t[0]);
            adj(1,3) = -t[1];
            adj(1,4) = s*(ad_bc*t[1] - (2*cd)*t[0]);
            adj(1,5) = s*((cc-dd)*t[0] - (ac-bd)*t[1]);

            zero(slice<2,0,4,2>(adj).instance());
        }
        
        template <class V>
        Vector<6,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<6>(v);
            const Matrix<2,2,Scalar> invA = latl::inverse(A);
            Matrix<2,2,Scalar> M;
            M(0,0) = v[3] + v[4];
            M(0,1) = v[5] - v[2];
            M(1,0) = v[5] + v[2];
            M(1,1) = v[3] - v[4];
            Matrix<2,2,Scalar> AMAinv = A * M * invA;
            
            Vector<6,Scalar> adj_v;
            slice<0,2>(adj_v) = A*slice<0,2>(v) - AMAinv*t;
            adj_v[2] = Scalar(0.5)*(AMAinv(1,0) - AMAinv(0,1));
            adj_v[3] = Scalar(0.5)*(AMAinv(0,0) + AMAinv(1,1));
            adj_v[4] = Scalar(0.5)*(AMAinv(0,0) - AMAinv(1,1));
            adj_v[5] = Scalar(0.5)*(AMAinv(1,0) + AMAinv(0,1));
            
            return adj_v;
        }

        template <class V>
        Vector<6,Scalar> adjoint_T_times(const AbstractVector<V>& v) const
        {
            assert_size_is<6>(v);
            Matrix<6,6,Scalar> adj;
            adjoint(adj);
            return adj.T() * v;
        }
    };
}

#endif
