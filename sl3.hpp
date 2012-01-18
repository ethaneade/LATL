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

#ifndef LATL_SL3_HPP
#define LATL_SL3_HPP

#include <latl/aff2.hpp>
#include <latl/lu.hpp>

namespace latl
{    
    template <class Scalar=double>
    class SL3 {
    private:
        Matrix<3,3,Scalar> H;
        explicit SL3(const Matrix<3,3,Scalar>& H_) : H(H_) {}

        template <class V>
        static void to_algebra(const AbstractVector<V>& v, Matrix<3,3,Scalar>& M) {
            assert_size_is<8>(v);
            M(0,2) = v[0];
            M(1,2) = v[1];
            M(0,0) = v[3] + v[4];
            M(0,1) = v[5] - v[2];
            M(1,0) = v[5] + v[2];
            M(1,1) = v[3] - v[4];
            M(2,0) = v[6];
            M(2,1) = v[7];
            M(2,2) = -2*v[3];
        }

        template <class V>
        static void from_algebra(const Matrix<3,3,Scalar>& M, AbstractVector<V>& v) {
            assert_size_is<8>(v);
            v[0] = M(0,2);
            v[1] = M(1,2);
            v[2] = Scalar(0.5)*(M(1,0) - M(0,1));
            v[3] = Scalar(0.5)*(M(0,0) + M(1,1));
            v[4] = Scalar(0.5)*(M(0,0) - M(1,1));
            v[5] = Scalar(0.5)*(M(1,0) + M(0,1));
            v[6] = M(2,0);
            v[7] = M(2,1);
        }        
        
        
    public:
        typedef Scalar scalar_type;
        enum {DoF = 8};
        
        SL3() : H(Identity()) {}

        SL3(const Uninitialized& u) {}

        template <class Mat>
        static SL3 from_matrix(const FixedMatrix<3,3,Mat>& H)
        {
            return SL3(H / latl::cbrt(determinant(H)));
        }
        
        SL3(const Aff2<Scalar>& aff2)
        {
            slice<0,0,2,2>(H) = aff2.linear() / latl::sqrt(determinant(aff2.linear()));
            H.T()[2] = unproject(aff2.translation());
            H(2,0) = H(2,1) = 0;
        }
        
        template <class T>
        explicit SL3(const SL3<T>& other) : H(other.H) {}

        void rectify() { H *= 1/latl::cbrt(determinant(H)); }
        
        const Matrix<3,3,Scalar>& matrix() const { return H; }
        
        template <class V>
        Vector<3,Scalar> operator*(const FixedVector<3,V>& x) const {
            return H*x;
        }
                
        SL3 operator*(const SL3& other) const {
            return SL3(H*other.H);
        }
        
        SL3 inverse() const {
            LU<3,Scalar> lu(H);
            Matrix<3,3,Scalar> Hinv;
            lu.get_inverse(Hinv);
            return SL3(Hinv);
        }

        template <class V>
        bool ln(FixedVector<8,V>& v) const {
            Matrix<3,3,Scalar> M;
            if (!log(H, M))
                return false;
            from_algebra(M, v);
            return true;
        }

        Vector<8,Scalar> ln() const {
            Vector<8,Scalar> v;
            bool success = ln(v);
            assert(success);
            return v;
        }

        template <class V>
        static SL3 exp(const AbstractVector<V>& v) {
            assert_size_is<8>(v);           
            Matrix<3,3,Scalar> M;
            to_algebra(v, M);
            Matrix<3,3,Scalar> H;
            latl::exp(M, H);            
            return SL3(H);
        }

        template <class Mat>
        void adjoint(AbstractMatrix<Mat>& adj) const
        {
            assert_shape_is<8,8>(adj);

            const Matrix<3,3,Scalar>& Hinv = inverse().H;
            
            Matrix<9,9,Scalar> C;
            slice<0,0,3,3>(C) = Hinv.T() * H(0,0);
            slice<0,3,3,3>(C) = Hinv.T() * H(0,1);
            slice<0,6,3,3>(C) = Hinv.T() * H(0,2);
            slice<3,0,3,3>(C) = Hinv.T() * H(1,0);
            slice<3,3,3,3>(C) = Hinv.T() * H(1,1);
            slice<3,6,3,3>(C) = Hinv.T() * H(1,2);
            slice<6,0,3,3>(C) = Hinv.T() * H(2,0);
            slice<6,3,3,3>(C) = Hinv.T() * H(2,1);
            slice<6,6,3,3>(C) = Hinv.T() * H(2,2);

            Matrix<9,8,Scalar> CE;
            CE.T()[0] = C.T()[2];
            CE.T()[1] = C.T()[5];
            CE.T()[2] = C.T()[3] - C.T()[1];
            CE.T()[3] = C.T()[0] + C.T()[4] - 2*C.T()[8];
            CE.T()[4] = C.T()[0] - C.T()[4];
            CE.T()[5] = C.T()[1] + C.T()[3];
            CE.T()[6] = C.T()[6];
            CE.T()[7] = C.T()[7];

            adj[0] = CE[2];
            adj[1] = CE[5];
            adj[2] = Scalar(0.5) * (CE[3] - CE[1]);
            adj[3] = Scalar(0.5) * (CE[0] + CE[4]);
            adj[4] = Scalar(0.5) * (CE[0] - CE[4]);
            adj[5] = Scalar(0.5) * (CE[3] + CE[1]);
            adj[6] = CE[6];
            adj[7] = CE[7];            
        }

        template <class V>
        Vector<8,Scalar> adjoint_times(const AbstractVector<V>& v) const {
            assert_size_is<8>(v);
            Matrix<3,3,Scalar> M;
            to_algebra(v, M);

            Vector<8,Scalar> adj_v;
            from_algebra(H * M * inverse().H, adj_v);
            return adj_v;
        }

        template <class V>
        Vector<8,Scalar> adjoint_T_times(const AbstractVector<V>& v) const {
            assert_size_is<8>(v);
            Vector<8,Scalar> sv = v;
            slice<2,4>(sv) *= Scalar(0.5);
            Matrix<3,3,Scalar> M;
            to_algebra(sv, M);
            M(2,2) = 0;

            Vector<8,Scalar> adj_v;
            M = H.T() * M * inverse().H.T();
            from_algebra(M, adj_v);
            adj_v[3] -= M(2,2);
            slice<2,4>(adj_v) *= 2;
            return adj_v;
        }

        // Jacobian of project(H * unproject(uv)) around H = I
        template <class V, class Mat>
        static void projection_jacobian(const AbstractVector<V>& uv,
                                        FixedMatrix<2,8,Mat>& J)
        {
            assert_size_is<2>(uv);
            Scalar u = uv[0], v = uv[1];
            J(0,0) = J(1,1) = 1;
            J(0,1) = J(1,0) = 0;
            J(0,2) = -v;
            J(1,2) = u;
            J(0,3) = 3*u;
            J(1,3) = 3*v;
            J(0,4) = u;
            J(1,4) = -v;
            J(0,5) = v;
            J(1,5) = u;
            J(0,6) = -u*u;
            J(1,6) = J(0,7) = -v*u;
            J(1,7) = -v*v;            
        }
                                        
    };
}

#endif
