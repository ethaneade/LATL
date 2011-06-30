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

#ifndef LATL_MATRIX_HPP
#define LATL_MATRIX_HPP

#include <latl/debug.hpp>
#include <latl/decl.hpp>
#include <latl/vector.hpp>

namespace latl
{

    //
    // AbstractMatrix definition
    //
    
    template <class Mat>
    struct AbstractMatrix {
	Mat& instance() { return static_cast<Mat&>(*this); }
	const Mat& instance() const { return static_cast<const Mat&>(*this); }

	typedef typename matrix_traits<Mat>::scalar_t scalar_t;
	typedef typename matrix_traits<Mat>::row_t row_t;
	typedef typename matrix_traits<Mat>::col_t col_t;
	typedef typename matrix_traits<Mat>::transpose_t transpose_t;

	enum {
            static_rows = matrix_traits<Mat>::static_rows,
            static_cols = matrix_traits<Mat>::static_cols,
            static_row_stride = matrix_traits<Mat>::static_row_stride,
            static_col_stride = matrix_traits<Mat>::static_col_stride,
        };

        int rows() const { return instance().mrows(); }
        int cols() const { return instance().mcols(); }
        
	const row_t operator[](int i) const {
            LATL_CHECK_BOUNDS(i, 0, rows());
            return instance().row(i);
        }
        
	row_t operator[](int i) {
            LATL_CHECK_BOUNDS(i, 0, rows());
            return instance().row(i);
        }

	scalar_t operator()(int i, int j) const {
            LATL_CHECK_BOUNDS(i, 0, rows());
            LATL_CHECK_BOUNDS(j, 0, cols());
            return instance().at(i,j);
        }
        
	scalar_t& operator()(int i, int j) {
            LATL_CHECK_BOUNDS(i, 0, rows());
            LATL_CHECK_BOUNDS(j, 0, cols());
            return instance().at(i,j);
        }
        
        transpose_t T() { return instance().transpose(); }
        const transpose_t T() const { return instance().transpose(); }
        
        scalar_t* data() { return instance().mdata(); }
        const scalar_t* data() const { return instance.mdata(); }

        template <class T>
        void assign(const AbstractMatrix<T>& other)
        {
            assert_same_shape(instance(), other.instance());
            unchecked_assign(other);
        }

        template <class T>
        void unchecked_assign(const AbstractMatrix<T>& other)
        {
            for (int i=0; i<rows(); ++i)
                for (int j=0; j<cols(); ++j)
                    (*this)(i,j) = other(i,j);
        }

        template <class E>
        Mat& operator=(const MatrixExpr<E>& e) {
            e.instance()(this->instance());
            return this->instance();
        }        
    };

    // 
    // Fill for constructors
    // 
    
    template <class Mat, class S>
    void fill(AbstractMatrix<Mat>& m, S s)
    {
        for (int i=0; i<m.rows(); ++i)
            for (int j=0; j<m.cols(); ++j)
                m(i,j) = s;
    }

    //
    // FixedMatrix definitions
    //
    
    template <int N, class Parent>
    struct FixedWidthMatrix : public Parent
    {
        int mcols() const { return N; }
    };

    template <int M, class Parent>
    struct FixedHeightMatrix : public Parent
    {
        int mrows() const { return M; }        
    };

    template <int M, int N, class Mat>
    struct FixedMatrix : public FixedHeightMatrix<M,
                                                  FixedWidthMatrix<N,
                                                                   AbstractMatrix<Mat> > >
    {
        int mrows() const { return M; }        
        int mcols() const { return N; }

    };

    template <class Mat>
    struct DynamicMatrix : public AbstractMatrix<Mat>
    {
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

    template <class Traits, int M = Traits::static_rows, int N = Traits::static_cols>
    struct MatrixBaseClass {
        typedef FixedMatrix<M, N, typename Traits::matrix_t> type;
    };

    template <class Traits, int N>
    struct MatrixBaseClass<Traits,-1,N> {
        typedef FixedWidthMatrix<N,AbstractMatrix<typename Traits::matrix_t> > type;
    };

    template <class Traits, int M>
    struct MatrixBaseClass<Traits,M,-1> {
        typedef FixedHeightMatrix<M,AbstractMatrix<typename Traits::matrix_t> > type;
    };

    template <class Traits>
    struct MatrixBaseClass<Traits,-1,-1> {
        typedef DynamicMatrix<typename Traits::matrix_t> type;
    };
    
    
    template <class Mat>
    struct TransposeMatrix : public MatrixBaseClass<matrix_traits<TransposeMatrix<Mat> > >::type
    {
    private:
        Mat *mat;
    public:
        typedef typename matrix_traits<TransposeMatrix<Mat> >::scalar_t scalar_t;
        typedef typename matrix_traits<TransposeMatrix<Mat> >::row_t row_t;
        typedef typename matrix_traits<TransposeMatrix<Mat> >::col_t col_t;

        TransposeMatrix(Mat& m)
        {
            mat = &m;
        }
        
        scalar_t* mdata() { return mat->mdata(); }
        const scalar_t* mdata() const { return mat->mdata(); }

        scalar_t at(int i, int j) const { return mat->at(j,i); }
        scalar_t& at(int i, int j) { return mat->at(j,i); }

        row_t row(int i) { return mat->col(i); }
        const row_t row(int i) const { return mat->col(i); }

        col_t col(int i) { return mat->row(i); }
        const col_t col(int i) const { return mat->row(i); }
        
        Mat& transpose() { return *mat; }
        const Mat& transpose() const { return *mat; }        
    };
        
    template <int M, int N, class Scalar>
    class Matrix : public FixedMatrix<M,N,Matrix<M,N,Scalar> >
    {
    private:
        Scalar x[M*N];
    public:
        typedef Matrix<M,N,Scalar> self_t;
        typedef matrix_traits<self_t> traits;
        typedef typename traits::row_t row_t;
        typedef typename traits::col_t col_t;
        typedef typename traits::transpose_t transpose_t;

        Matrix() {}
        
        Matrix(Scalar s) {
            fill(*this, s);
        }
        
        Scalar* mdata() { return x; }
        const Scalar* mdata() const { return x; }
        
        Scalar at(int i, int j) const { return x[i*N + j]; }
        Scalar& at(int i, int j) { return x[i*N + j]; }

        row_t row(int i) { return row_t(x+i*N); }
        const row_t row(int i) const { return row_t(const_cast<Scalar*>(x+i*N)); }
        
        col_t col(int i) { return col_t(x+i); }
        const col_t col(int i) const { return col_t(const_cast<Scalar*>(x+i)); }

        transpose_t transpose() {
            return transpose_t(*this);
        }
        
        const transpose_t transpose() const {
            return transpose_t(const_cast<self_t&>(*this));
        }
    };    
}

#endif
