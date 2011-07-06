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

#ifndef LATL_MATRIX_OPS_HPP
#define LATL_MATRIX_OPS_HPP

#include <latl/debug.hpp>
#include <latl/matrix.hpp>
#include <latl/expr.hpp>
#include <latl/scalar.hpp>
#include <latl/vector_ops.hpp>

namespace latl
{
#define LATL_MS(M) typename matrix_traits<M>::scalar_t
#define LATL_WIDER(A,B) typename Wider<A,B>::type
#define LATL_WIDER_MS(M1,M2) typename Wider<LATL_MS(M1),LATL_MS(M2)>::type
    
    // Scalar * Matrix
    
    template <class S, class Mat>
    struct ScalarMatrixProduct : public MatrixExpr<ScalarMatrixProduct<S,Mat> > {
        const AbstractMatrix<Mat>& m;
        S s;

        ScalarMatrixProduct(const AbstractMatrix<Mat>& m_, S s_) : m(m_), s(s_) {}
        int rows() const { return m.rows(); }
        int cols() const { return m.cols(); }
        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            assert_same_shape(q, m);
            for (int i=0; i<m.rows(); ++i)
                for (int j=0; j<m.cols(); ++j)
                    q(i,j) = m(i,j) * s;
        }

        typedef Matrix<matrix_traits<Mat>::static_rows,
                       matrix_traits<Mat>::static_cols,
                       typename Wider<S,LATL_MS(Mat)>::type> result_t;
    };
   
    template <class S, class Mat>
    typename ScalarMatrixProduct<typename ScalarType<S>::type, Mat>::result_t
    operator*(S s, const AbstractMatrix<Mat>& m) {
        return ScalarMatrixProduct<S,Mat>(m,s);
    }

    template <class S, class Mat>
    typename ScalarMatrixProduct<typename ScalarType<S>::type, Mat>::result_t
    operator*(const AbstractMatrix<Mat>& m, S s) {
        return ScalarMatrixProduct<S,Mat>(m,s);
    }

    //
    // Matrix / Scalar 
    //

    
    // For invertible types, multiply by the reciprocal
    template <class S, class Mat>    
    typename ScalarMatrixProduct<typename ScalarType<LATL_WIDER(S,LATL_MS(Mat))>::recip_type,Mat>::result_t
    operator/(const AbstractMatrix<Mat>& m, S s) {
        typedef LATL_WIDER(S,LATL_MS(Mat)) W;
        return ScalarMatrixProduct<W,Mat>(m,1/W(s));
    }    

    template <class Mat, class S>
    struct MatrixQuotient : public MatrixExpr<MatrixQuotient<Mat,S> > {
        const AbstractMatrix<Mat>& m;
        S s;

        MatrixQuotient(const AbstractMatrix<Mat>& m_, S s_) : m(m_), s(s_) {}
        int rows() const { return m.rows(); }
        int cols() const { return m.cols(); }
        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            assert_same_shape(q, m);
            for (int i=0; i<m.rows(); ++i)
                for (int j=0; j<m.cols(); ++j)
                    q(i,j) = m(i,j)/s;
        }

        typedef Matrix<matrix_traits<Mat>::static_rows,
                       matrix_traits<Mat>::static_cols,
                       typename Wider<S,LATL_MS(Mat)>::type> result_t;
    };

    // For non-invertible types, divide each element
    template <class S, class Mat>
    typename MatrixQuotient<Mat,typename ScalarType<LATL_WIDER(S,LATL_MS(Mat))>::nonrecip_type>::result_t
    operator/(const AbstractMatrix<Mat>& m, S s) {
        typedef LATL_WIDER(S,LATL_MS(Mat)) W;
        return MatrixQuotient<Mat,W>(m,s);
    }

    
    //
    // Matrix negation
    //
     
    template <class Mat>
    struct MatrixNegation : public MatrixExpr<MatrixNegation<Mat> > {
        const AbstractMatrix<Mat>& m;
        int rows() const { return m.rows(); }
        int cols() const { return m.cols(); }
        MatrixNegation(const AbstractMatrix<Mat>& m_) : m(m_) {}
        
        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            assert_same_shape(q, m);
            for (int i=0; i<m.rows(); ++i)
                for (int j=0; j<m.cols(); ++j)
                    q(i,j) = -m(i,j);
        }
        typedef Matrix<matrix_traits<Mat>::static_rows,
                       matrix_traits<Mat>::static_cols,
                       LATL_MS(Mat)> result_t;
    };
    
    template <class Mat>
    typename MatrixNegation<Mat>::result_t
    operator-(const AbstractMatrix<Mat>& m) {
        return MatrixNegation<Mat>(m);
    }

    //
    // Matrix sum
    //
     
    template <class Mat1, class Mat2>
    struct MatrixSum : public MatrixExpr<MatrixSum<Mat1,Mat2> > {
        const AbstractMatrix<Mat1>& a;
        const AbstractMatrix<Mat2>& b;
        int rows() const { return a.rows(); }
        int cols() const { return a.cols(); }
        MatrixSum(const AbstractMatrix<Mat1>& a_,
                  const AbstractMatrix<Mat2>& b_)
            : a(a_), b(b_)
        {
            assert_same_shape(a,b);
        }
        
        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            assert_same_shape(q, a);
            for (int i=0; i<a.rows(); ++i)
                for (int j=0; j<a.cols(); ++j)
                    q(i,j) = a(i,j) + b(i,j);
        }

        typedef Matrix<MaxInt<matrix_traits<Mat1>::static_rows, matrix_traits<Mat2>::static_rows>::value,
                       MaxInt<matrix_traits<Mat1>::static_cols, matrix_traits<Mat2>::static_cols>::value,
                       LATL_WIDER_MS(Mat1,Mat2)> result_t;
    };
    
    template <class Mat1, class Mat2>
    typename MatrixSum<Mat1,Mat2>::result_t
    operator+(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatrixSum<Mat1,Mat2>(a,b);
    }

    //
    // Matrix difference
    //
     
    template <class Mat1, class Mat2>
    struct MatrixDiff : public MatrixExpr<MatrixDiff<Mat1,Mat2> > {
        const AbstractMatrix<Mat1>& a;
        const AbstractMatrix<Mat2>& b;
        int rows() const { return a.rows(); }
        int cols() const { return a.cols(); }
        MatrixDiff(const AbstractMatrix<Mat1>& a_,
                   const AbstractMatrix<Mat2>& b_)
            : a(a_), b(b_)
        {
            assert_same_shape(a,b);
        }
        
        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            assert_same_shape(q, a);
            for (int i=0; i<a.rows(); ++i)
                for (int j=0; j<a.cols(); ++j)
                    q(i,j) = a(i,j) - b(i,j);
        }

        typedef Matrix<MaxInt<matrix_traits<Mat1>::static_rows, matrix_traits<Mat2>::static_rows>::value,
                       MaxInt<matrix_traits<Mat1>::static_cols, matrix_traits<Mat2>::static_cols>::value,
                       LATL_WIDER_MS(Mat1,Mat2)> result_t;
    };
    
    template <class Mat1, class Mat2>
    typename MatrixDiff<Mat1,Mat2>::result_t
    operator-(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatrixDiff<Mat1,Mat2>(a,b);
    }


    //
    // Matrix * Vector
    //

    template <class Mat, class V>
    struct MatVecProduct  : public VectorExpr<MatVecProduct<Mat,V> > {
        const AbstractMatrix<Mat>& m;
        const AbstractVector<V>& v;

        int size() const { return m.rows(); }
        
        MatVecProduct(const AbstractMatrix<Mat>& m_,
                      const AbstractVector<V>& v_)
            : m(m_), v(v_)
        {
            CheckEquality<matrix_traits<Mat>::static_cols,
                vector_traits<V>::static_size>::eval(m.cols(), v.size());
        }
        
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            CheckEquality<matrix_traits<Mat>::static_rows,
                vector_traits<W>::static_size>::eval(m.rows(), v.size());
        
            for (int i=0; i<m.rows(); ++i)
                w[i] = m[i] * v;
        }

        typedef Vector<matrix_traits<Mat>::static_rows,
                       typename Wider<typename matrix_traits<Mat>::scalar_t,
                                      typename vector_traits<V>::scalar_t>::type> result_t;
    };

    template <class Mat, class V>
    typename MatVecProduct<Mat,V>::result_t
    operator*(const AbstractMatrix<Mat>& m, const AbstractVector<V>& v) {
        return MatVecProduct<Mat,V>(m,v);
    }
    
}

#endif
