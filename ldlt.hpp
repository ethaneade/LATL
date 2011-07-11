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

#ifndef LATL_LDLT_HPP
#define LATL_LDLT_HPP

#include <latl/ldlt.hpp>

namespace latl
{
    template <int N, class Scalar>
    class LDLT
    {
    public:
        MatrixHolder<N,N,Scalar> L;
        VectorHolder<N,Scalar> inv_diag;
        bool full_rank;

    public:
        LDLT() {
            full_rank = false;
        }
        
        template <class Mat>
        LDLT(const AbstractMatrix<Mat>& m) { compute(m); }
        

        template <class Mat>
        bool compute(const AbstractMatrix<Mat>& m) {
            assert_square(m);
            int n = m.rows();
            L.resize(n, n);
            assert_same_shape(L(), m);
            
            inv_diag.resize(n);
            full_rank = false;
            
            for (int i=0; i<n; ++i) {
                Scalar d = m(i,i);
                for (int k=0; k<i; ++k)
                    d -= L()(i,k) * L()(k,i);
                if (d == 0)
                    return false;
                L()(i,i) = d;
                inv_diag()[i] = 1/d;
                
                for (int j=i+1; j<n; ++j) {
                    Scalar y = m(i,j);
                    for (int k=0; k<i; ++k)
                        y -= L()(i,k) * L()(k,j);
                    L()(j,i) = inv_diag()[i] * y;
                    L()(i,j) = y;
                }
            }
            full_rank = true;
            return true;
        }

        template <class V>
        typename VectorHolder<N,Scalar>::type
        inverse_times(const AbstractVector<V>& v) const {
            assert_same_size(inv_diag(), v);
            const int n = L().rows();
            VectorHolder<N,Scalar> y(n);
            for (int i=0; i<n; ++i) {
                Scalar yi = v[i];
                for (int j=0; j<i; ++j)
                    yi -= L()(i, j) * y()[j];
                y()[i] = yi;
            }
            for (int i=n-1; i>=0; --i) {
                Scalar yi = inv_diag()[i] * y()[i];
                for (int j=i+1; j<n; ++j)
                    yi -= L()(j,i) * y()[j];
                y()[i] = yi;
            }
            return y();
        }

        bool is_full_rank() const { return full_rank; }
        
        Scalar determinant() const {
            if (!is_full_rank())
                return 0;
            Scalar det = L()(0,0);
            for (int i=1; i<L().rows(); ++i)
                det *= L()(i,i);
            return det;
        }

        template <class Mat>
        void get_Cholesky(AbstractMatrix<Mat>& chol_L)
        {
            assert(is_full_rank());
            assert_same_shape(L(), chol_L);
            
            const int n = L().rows();
            for (int j=0; j<n; ++j) {
                assert(L()(j,j) >= 0);
                Scalar sd = latl::sqrt(L()(j,j));
                chol_L(j,j) = sd;
                for (int i=j+1; i<n; ++i) {
                    chol_L(i,j) = sd * L()(i,j);
                    chol_L(j,i) = 0;
                }
            }
        }
    };
}

#endif
