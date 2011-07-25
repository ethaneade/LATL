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

#ifndef LATL_MATRIX_FUNCS_HPP
#define LATL_MATRIX_FUNCS_HPP

#include <latl/matrix.hpp>
#include <latl/scalar.hpp>

namespace latl
{

    template <class Mat>
    void transpose(AbstractMatrix<Mat>& m)
    {
        assert_square(m);
        for (int i=0; i<m.rows(); ++i) {
            for (int j=i+1; j<m.cols(); ++j) {
                LATL_MS(Mat) tmp = m(i,j);
                m(i,j) = m(j,i);
                m(j,i) = tmp;
            }
        }
    }
    
    template <class Mat>
    LATL_MS(Mat) trace(const AbstractMatrix<Mat>& m)
    {
        assert_square(m);
        LATL_MS(Mat) sum = 0;
        for (int i=0; i<m.rows(); ++i)
            sum += m(i,i);
        return sum;
    }
    
    
    template <class Mat>
    void orthonormalize(AbstractMatrix<Mat>& m) {
        for (int i=0; i<m.rows(); ++i) {
            normalize(m[i].instance());
            for (int j=i+1; j<m.rows(); ++j)
                m[j] -= (m[j] * m[i]) * m[i];
        }
    }

    template <class Mat>
    Matrix<2,2,LATL_MS(Mat)> inverse(const FixedMatrix<2,2,Mat>& m)
    {
        typedef LATL_MS(Mat) Scalar;
        Scalar det = m(0,0)*m(1,1) - m(0,1)*m(1,0);
        assert(det != 0);

        Scalar inv_det = 1 / det;
        Matrix<2,2,Scalar> inv;
        inv(0,0) = inv_det * m(1,1);
        inv(1,1) = inv_det * m(0,0);
        inv(0,1) = inv_det * -m(0,1);
        inv(1,0) = inv_det * -m(1,0);
        return inv;
    }

    template <class V>
    Matrix<3,3,LATL_VS(V)> skew_matrix(const AbstractVector<V>& w)
    {
        assert_size_is<3>(w);
        Matrix<3,3,LATL_VS(V)> s;
        s(0,0) = s(1,1) = s(2,2) = 0;
        s(0,1) = -w[2];
        s(0,2) = w[1];
        s(1,0) = w[2];
        s(1,2) = -w[0];
        s(2,0) = -w[1];
        s(2,1) = w[0];
        return s;
    }

    template <class V, class Mat, class Op>
    void outer_product_upper(const AbstractVector<V>& v, const Op& op, AbstractMatrix<Mat>& m)
    {
        assert_same_size(v, m[0]);
        assert_same_size(v, m.T()[0]);
        typedef LATL_VS(V) Scalar;
        const int n = v.size();
        for (int i=0; i<n; ++i) {
            Scalar vi = v[i];
            for (int j=i; j<n; ++j)
                op(m(i,j), vi*v[j]);
        }
    }

    template <class V, class Mat, class Op>
    void outer_product_upper(const AbstractVector<V>& v, LATL_VS(V) s, const Op& op, AbstractMatrix<Mat>& m)
    {
        assert_same_size(v, m[0]);
        assert_same_size(v, m.T()[0]);
        typedef LATL_VS(V) Scalar;
        const int n = v.size();
        for (int i=0; i<n; ++i) {
            Scalar svi = s*v[i];
            for (int j=i; j<n; ++j)
                op(m(i,j), svi*v[j]);
        }
    }
    
}

#endif
