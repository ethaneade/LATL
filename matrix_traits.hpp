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

#ifndef LATL_MATRIX_TRAITS_HPP
#define LATL_MATRIX_TRAITS_HPP

#include <latl/decl.hpp>

namespace latl {

    template <class Mat>
    struct matrix_traits<AbstractMatrix<Mat> > {
        typedef matrix_traits<Mat> traits_t;
        typedef Mat matrix_t;
        typedef typename traits_t::scalar_t scalar_t;
        typedef typename traits_t::row_t row_t;
        typedef typename traits_t::col_t col_t;
        typedef typename traits_t::transpose_t transpose_t;

        enum {
            static_rows = traits_t::static_rows,
            static_cols = traits_t::static_cols,
            static_row_stride = traits_t::static_row_stride,
            static_col_stride = traits_t::static_col_stride,
        };
    };
    
    template <int M, int N, class Scalar>
    struct matrix_traits<Matrix<M,N,Scalar> > {
        typedef Matrix<M,N,Scalar> matrix_t;
        typedef Scalar scalar_t;
        typedef RefVector<N,1,Scalar> row_t;
        typedef RefVector<M,N,Scalar> col_t;
        typedef TransposeMatrix<matrix_t> transpose_t;

        enum {
            static_rows = M,
            static_cols = N,
            static_row_stride = N,
            static_col_stride = 1,
        };
    };
    
    template <class Mat>
    struct matrix_traits<TransposeMatrix<Mat> > {
        typedef TransposeMatrix<Mat> matrix_t;
        typedef typename matrix_traits<Mat>::scalar_t scalar_t;
        typedef typename matrix_traits<Mat>::col_t row_t;
        typedef typename matrix_traits<Mat>::row_t col_t;
        typedef Mat& transpose_t;

        enum {
            static_rows = matrix_traits<Mat>::static_cols,
            static_cols = matrix_traits<Mat>::static_rows,
            static_row_stride = matrix_traits<Mat>::static_col_stride,
            static_col_stride = matrix_traits<Mat>::static_row_stride,
        };
    };

    template <int M, int N, int RS, int CS, class Scalar>
    struct matrix_traits<RefMatrix<M,N,RS,CS,Scalar> > {
        typedef RefMatrix<M,N,RS,CS,Scalar> matrix_t;
        typedef Scalar scalar_t;
        typedef RefVector<N,CS,Scalar> row_t;
        typedef RefVector<M,RS,Scalar> col_t;
        typedef TransposeMatrix<matrix_t> transpose_t;

        enum {
            static_rows = M,
            static_cols = N,
            static_row_stride = RS,
            static_col_stride = CS,
        };
    };
    
}

#endif
