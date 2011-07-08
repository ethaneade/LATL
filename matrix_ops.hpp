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

    template <class Mat> template <class S>
    Mat& AbstractMatrix<Mat>::operator*=(S s) {
        for (int i=0; i<rows(); ++i)
            for (int j=0; j<cols(); ++j)
                (*this)(i,j) *= s;
        return instance();
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

    template <class W, class Mat>
    void divide(AbstractMatrix<Mat>& m, typename ScalarType<W>::recip_type s) {
        W r = 1/s;
        for (int i=0; i<m.rows(); ++i)
            for (int j=0; j<m.cols(); ++j)
                m(i,j) *= r;
    }

    template <class W, class Mat>
    void divide(AbstractMatrix<Mat>& m, typename ScalarType<W>::nonrecip_type s) {
        for (int i=0; i<m.rows(); ++i)
            for (int j=0; j<m.cols(); ++j)
                m(i,j) /= s;
    }
    
    template <class Mat> template <class S>
    Mat& AbstractMatrix<Mat>::operator/=(S s) {
        typedef typename Wider<LATL_MS(Mat),S>::type W;
        divide<W>(*this, s);
        return instance();
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

    template <class Mat1, class Mat2>
    struct BinOpResult {
        enum {
            M1 = matrix_traits<Mat1>::static_rows,
            M2 = matrix_traits<Mat2>::static_rows,
            N1 = matrix_traits<Mat1>::static_cols,
            N2 = matrix_traits<Mat2>::static_cols,
            M = MaxInt<M1,M2>::value,
            N = MaxInt<N1,N2>::value,
        };

        typedef Matrix<M,N, LATL_WIDER_MS(Mat1,Mat2)> type;
    };
    
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

        typedef typename BinOpResult<Mat1,Mat2>::type result_t;
    };
    
    template <class Mat1, class Mat2>
    typename MatrixSum<Mat1,Mat2>::result_t
    operator+(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatrixSum<Mat1,Mat2>(a,b);
    }

    template <class M1> template <class M2>
    M1& AbstractMatrix<M1>::operator+=(const AbstractMatrix<M2>& m) {
        assert_same_shape(*this, m);
        for (int i=0; i<rows(); ++i)
            for (int j=0; j<cols(); ++j)
                (*this)(i,j) += m(i,j);
        return instance();
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
        
        typedef typename BinOpResult<Mat1,Mat2>::type result_t;
    };
    
    template <class Mat1, class Mat2>
    typename MatrixDiff<Mat1,Mat2>::result_t
    operator-(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatrixDiff<Mat1,Mat2>(a,b);
    }

    template <class M1> template <class M2>
    M1& AbstractMatrix<M1>::operator-=(const AbstractMatrix<M2>& m) {
        assert_same_shape(*this, m);
        for (int i=0; i<rows(); ++i)
            for (int j=0; j<cols(); ++j)
                (*this)(i,j) -= m(i,j);
        return instance();
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

    //
    // Matrix * Matrix
    //

    template <class M1, class M2>
    struct MatrixProduct  : public MatrixExpr<MatrixProduct<M1,M2> > {
        const AbstractMatrix<M1>& a;
        const AbstractMatrix<M2>& b;

        int rows() const { return a.rows(); }
        int cols() const { return b.cols(); }
        
        MatrixProduct(const AbstractMatrix<M1>& a_,
                      const AbstractMatrix<M2>& b_)
            : a(a_), b(b_)
        {
            CheckEquality<matrix_traits<M1>::static_cols,
                matrix_traits<M2>::static_rows>::eval(a.cols(), b.rows());
        }
        
        template <class T>
        void operator()(AbstractMatrix<T>& c) const {
            CheckEquality<matrix_traits<T>::static_rows,
                matrix_traits<M1>::static_rows>::eval(c.rows(), a.rows());
            CheckEquality<matrix_traits<T>::static_cols,
                matrix_traits<M2>::static_cols>::eval(c.cols(), b.cols());
        
            for (int i=0; i<b.cols(); ++i) {
                c.col(i) = a * b.col(i);
            }
        }

        typedef typename Wider<typename matrix_traits<M1>::scalar_t,
                               typename matrix_traits<M2>::scalar_t>::type wider_t;

        typedef Matrix<matrix_traits<M1>::static_rows,
                       matrix_traits<M2>::static_cols,
                       wider_t> result_t;
    };

    template <class M1, class M2>
    typename MatrixProduct<M1,M2>::result_t
    operator*(const AbstractMatrix<M1>& a, const AbstractMatrix<M2>& b) {
        return MatrixProduct<M1,M2>(a,b);
    }

    //
    // Outer product
    //

    template <class V1, class V2>
    struct OuterProduct  : public MatrixExpr<OuterProduct<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;

        int rows() const { return a.size(); }
        int cols() const { return b.size(); }
        
        OuterProduct(const AbstractVector<V1>& a_,
                      const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
        }
        
        template <class T>
        void operator()(AbstractMatrix<T>& c) const {
            CheckEquality<matrix_traits<T>::static_rows,
                vector_traits<V1>::static_size>::eval(c.rows(), a.size());
            CheckEquality<matrix_traits<T>::static_cols,
                vector_traits<V2>::static_size>::eval(c.cols(), b.size());
        
            for (int i=0; i<a.size(); ++i)
                c[i] = a[i] * b;
        }

        typedef typename Wider<typename vector_traits<V1>::scalar_t,
                               typename vector_traits<V2>::scalar_t>::type wider_t;

        typedef Matrix<vector_traits<V1>::static_size,
                       vector_traits<V2>::static_size,
                       wider_t> result_t;
    };
    
    template <class V1, class V2>
    typename OuterProduct<V1,V2>::result_t
    outer_product(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return OuterProduct<V1,V2>(a,b);
    }


    //
    // (Vector | Vector) columnar stacking
    //

    template <class V1, class V2>
    struct VectorColCat : public MatrixExpr<VectorColCat<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = MaxInt<N1,N2>::value
        };
        int rows() const { return a.size(); }
        int cols() const { return 2; }
        VectorColCat(const AbstractVector<V1>& a_,
                     const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
            assert_same_size(a,b);
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_rows, N>::eval(q.rows(), rows());
            CheckEquality<matrix_traits<T>::static_cols, 2>::eval(q.cols(), 2);

            q.col(0) = a;
            q.col(1) = b;
        }

        typedef Matrix<N, 2, LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorColCat<V1,V2>::result_t
    operator|(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorColCat<V1,V2>(a,b);
    }

    //
    // [Vector]
    // [Vector] row-wise stacking
    //

    template <class V1, class V2>
    struct VectorRowCat : public MatrixExpr<VectorRowCat<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = MaxInt<N1,N2>::value
        };
        int rows() const { return 2; }
        int cols() const { return a.size(); }
        
        VectorRowCat(const AbstractVector<V1>& a_,
                     const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
            assert_same_size(a,b);
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_cols, N>::eval(q.cols(), cols());
            CheckEquality<matrix_traits<T>::static_rows, 2>::eval(q.rows(), 2);

            q[0] = a;
            q[1] = b;
        }

        typedef Matrix<2, N, LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorRowCat<V1,V2>::result_t
    operator%(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorRowCat<V1,V2>(a,b);
    }

    //
    // (Matrix | Vector) columnar append
    //

    template <class Mat, class V>
    struct MatVecColCat : public MatrixExpr<MatVecColCat<Mat,V> > {
        const AbstractMatrix<Mat>& m;
        const AbstractVector<V>& v;
        enum {
            M1 = matrix_traits<Mat>::static_rows,
            M2 = vector_traits<V>::static_size,
            M = MaxInt<M1,M2>::value,
            N1 = matrix_traits<Mat>::static_cols,
            N = (N1 == -1 ? -1 : N1 + 1),
        };
        int rows() const { return m.rows(); }
        int cols() const { return m.cols() + 1; }
        MatVecColCat(const AbstractMatrix<Mat>& m_,
                     const AbstractVector<V>& v_)
            : m(m_), v(v_)
        {
            CheckEquality<M1,M2>::eval(m.rows(), v.size());
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_rows, M>::eval(q.rows(), rows());
            CheckEquality<matrix_traits<T>::static_cols, N>::eval(q.cols(), cols());

            for (int i=0; i<m.cols(); ++i)
                q.col(i) = m.col(i);
            q.col(m.cols()) = v;
        }

        typedef Matrix<M, N, LATL_WIDER(LATL_MS(Mat),LATL_VS(V))> result_t;
    };

    template <class Mat, class V>
    typename MatVecColCat<Mat,V>::result_t
    operator|(const AbstractMatrix<Mat>& m, const AbstractVector<V>& v) {
        return MatVecColCat<Mat,V>(m,v);
    }

    //
    // [Matrix]
    // [Vector] row-wise append
    //

    template <class Mat, class V>
    struct MatVecRowCat : public MatrixExpr<MatVecRowCat<Mat,V> > {
        const AbstractMatrix<Mat>& m;
        const AbstractVector<V>& v;
        enum {
            M1 = matrix_traits<Mat>::static_rows,
            M = (M1 == -1 ? -1 : M1 + 1),
            N1 = matrix_traits<Mat>::static_cols,
            N2 = vector_traits<V>::static_size,
            N = MaxInt<N1,N2>::value,
        };
        int rows() const { return m.rows() + 1; }
        int cols() const { return m.cols(); }
        MatVecRowCat(const AbstractMatrix<Mat>& m_,
                     const AbstractVector<V>& v_)
            : m(m_), v(v_)
        {
            CheckEquality<N1,N2>::eval(m.cols(), v.size());
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_rows, M>::eval(q.rows(), rows());
            CheckEquality<matrix_traits<T>::static_cols, N>::eval(q.cols(), cols());

            for (int i=0; i<m.rows(); ++i)
                q[i] = m[i];
            q[m.rows()] = v;
        }

        typedef Matrix<M, N, LATL_WIDER(LATL_MS(Mat),LATL_VS(V))> result_t;
    };

    template <class Mat, class V>
    typename MatVecRowCat<Mat,V>::result_t
    operator%(const AbstractMatrix<Mat>& m, const AbstractVector<V>& v) {
        return MatVecRowCat<Mat,V>(m,v);
    }

    //
    // (Matrix | Matrix) side-by-side append
    //

    template <class Mat1, class Mat2>
    struct MatMatColCat : public MatrixExpr<MatMatColCat<Mat1,Mat2> > {
        const AbstractMatrix<Mat1>& a;
        const AbstractMatrix<Mat2>& b;
        enum {
            M1 = matrix_traits<Mat1>::static_rows,
            M2 = matrix_traits<Mat2>::static_rows,
            M = MaxInt<M1,M2>::value,
            N1 = matrix_traits<Mat1>::static_cols,
            N2 = matrix_traits<Mat2>::static_cols,
            N = (N1 == -1 || N2 == -1) ? -1 : N1 + N2,
        };
        int rows() const { return a.rows(); }
        int cols() const { return a.cols() + b.cols(); }
        MatMatColCat(const AbstractMatrix<Mat1>& a_,
                     const AbstractMatrix<Mat2>& b_)
            : a(a_), b(b_)
        {
            CheckEquality<M1,M2>::eval(a.rows(), b.rows());
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_rows, M>::eval(q.rows(), rows());
            CheckEquality<matrix_traits<T>::static_cols, N>::eval(q.cols(), cols());

            for (int i=0; i<a.cols(); ++i)
                q.col(i) = a.col(i);
            for (int i=0; i<b.cols(); ++i)
                q.col(a.cols() + i) = b.col(i);
        }

        typedef Matrix<M, N, LATL_WIDER_MS(Mat1,Mat2)> result_t;
    };

    template <class Mat1, class Mat2>
    typename MatMatColCat<Mat1,Mat2>::result_t
    operator|(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatMatColCat<Mat1,Mat2>(a,b);
    }

    //
    // [Matrix]
    // [Matrix] stacking
    //

    template <class Mat1, class Mat2>
    struct MatMatRowCat : public MatrixExpr<MatMatRowCat<Mat1,Mat2> > {
        const AbstractMatrix<Mat1>& a;
        const AbstractMatrix<Mat2>& b;
        enum {
            M1 = matrix_traits<Mat1>::static_rows,
            M2 = matrix_traits<Mat2>::static_rows,
            N1 = matrix_traits<Mat1>::static_cols,
            N2 = matrix_traits<Mat2>::static_cols,
            M = (M1 == -1 || M2 == -1) ? -1 : M1 + M2,
            N = MaxInt<N1,N2>::value,
        };
        int rows() const { return a.rows(); }
        int cols() const { return a.cols() + b.cols(); }
        MatMatRowCat(const AbstractMatrix<Mat1>& a_,
                     const AbstractMatrix<Mat2>& b_)
            : a(a_), b(b_)
        {
            CheckEquality<M1,M2>::eval(a.rows(), b.rows());
        }

        template <class T>
        void operator()(AbstractMatrix<T>& q) const {
            CheckEquality<matrix_traits<T>::static_rows, M>::eval(q.rows(), rows());
            CheckEquality<matrix_traits<T>::static_cols, N>::eval(q.cols(), cols());

            for (int i=0; i<a.rows(); ++i)
                q[i] = a[i];
            for (int i=0; i<b.rows(); ++i)
                q[a.rows() + i] = b[i];
        }

        typedef Matrix<M, N, LATL_WIDER_MS(Mat1,Mat2)> result_t;
    };

    template <class Mat1, class Mat2>
    typename MatMatRowCat<Mat1,Mat2>::result_t
    operator%(const AbstractMatrix<Mat1>& a, const AbstractMatrix<Mat2>& b) {
        return MatMatRowCat<Mat1,Mat2>(a,b);
    }
    

    
}

#endif
