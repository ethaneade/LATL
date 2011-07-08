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
#include <latl/vector.hpp>
#include <latl/matrix_traits.hpp>

namespace latl
{    
    template <class A, class B>
    void assert_same_shape(const AbstractMatrix<A>& a, const AbstractMatrix<B>& b) {
        CheckEquality<matrix_traits<A>::static_rows,
            matrix_traits<B>::static_rows>::eval(a.rows(), b.rows());
        CheckEquality<matrix_traits<A>::static_cols,
            matrix_traits<B>::static_cols>::eval(a.cols(), b.cols());
    }

    template <class A>
    void assert_square(const AbstractMatrix<A>& a) {
        CheckEquality<matrix_traits<A>::static_rows,
            matrix_traits<A>::static_cols>::eval(a.rows(), a.cols());
    }
    
    
    //
    // AbstractMatrix definition
    //
    
    template <class Mat>
    class AbstractMatrix {
    public:
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
        int row_stride() const { return instance().mrow_stride(); }
        int col_stride() const { return instance().mcol_stride(); }        
        
	const row_t operator[](int i) const {
            LATL_CHECK_BOUNDS(i, 0, rows());
            return instance().row(i);
        }
        
	row_t operator[](int i) {
            LATL_CHECK_BOUNDS(i, 0, rows());
            return instance().row(i);
        }

	const col_t col(int i) const {
            LATL_CHECK_BOUNDS(i, 0, cols());
            return instance().col(i);
        }
        
	col_t col(int i) {
            LATL_CHECK_BOUNDS(i, 0, cols());
            return instance().col(i);
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
        typename Constify<transpose_t>::type T() const { return instance().transpose(); }
        
        scalar_t* data() { return instance().mdata(); }
        const scalar_t* data() const { return instance.mdata(); }

        scalar_t* data_at(int i, int j) {
            return data() + row_stride()*i + col_stride()*j;
        }

        const scalar_t* data_at(int i, int j) const {
            return data() + row_stride()*i + col_stride()*j;
        }
        
        template <class T>
        void assign(const AbstractMatrix<T>& other)
        {
            assert_same_shape(instance(), other.instance());
            unchecked_assign(other);
        }

        template <class E>
        Mat& operator=(const MatrixExpr<E>& e) {
            e.instance()(this->instance());
            return this->instance();
        }

        template <class T> Mat& operator+=(const AbstractMatrix<T>& m);
        template <class T> Mat& operator-=(const AbstractMatrix<T>& m);
        template <class S> Mat& operator*=(S s);
        template <class S> Mat& operator/=(S s);
        
    private:
        template <class T>
        void unchecked_assign(const AbstractMatrix<T>& other)
        {
            for (int i=0; i<rows(); ++i)
                for (int j=0; j<cols(); ++j)
                    (*this)(i,j) = other(i,j);
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
        typedef TransposeMatrix<Mat> self_t;
        typedef typename matrix_traits<self_t>::scalar_t scalar_t;
        typedef typename matrix_traits<self_t>::row_t row_t;
        typedef typename matrix_traits<self_t>::col_t col_t;

        TransposeMatrix(Mat& m)
        {
            mat = &m;
        }

        template <class T>
        self_t& operator=(const AbstractMatrix<T>& other) {
            assign(other);
            return *this;
        }
        
        scalar_t* mdata() { return mat->mdata(); }
        const scalar_t* mdata() const { return mat->mdata(); }

        int mrows() const { return mat->mcols(); }
        int mcols() const { return mat->mrows(); }
        int mrow_stride() const { return mat->mcol_stride(); }
        int mcol_stride() const { return mat->mrow_stride(); }

        scalar_t at(int i, int j) const { return mat->at(j,i); }
        scalar_t& at(int i, int j) { return mat->at(j,i); }

        row_t row(int i) { return mat->col(i); }
        const row_t row(int i) const { return mat->col(i); }

        col_t col(int i) { return mat->row(i); }
        const col_t col(int i) const { return mat->row(i); }
        
        Mat& transpose() { return *mat; }
        const Mat& transpose() const { return *mat; }        
    };
        
#define LATL_REFM_COMMON()                      \
        typedef typename matrix_traits<self_t>::row_t row_t;            \
        typedef typename matrix_traits<self_t>::col_t col_t;            \
        typedef typename matrix_traits<self_t>::transpose_t transpose_t; \
                                                                        \
        RefMatrix(Scalar *x_) : x(x_) {}                                \
        RefMatrix& operator=(const RefMatrix& other) {                  \
            assign(other);                                              \
            return *this;                                               \
        }                                                               \
        template <class T>                                              \
        RefMatrix& operator=(const AbstractMatrix<T>& other) {          \
            assign(other);                                              \
            return *this;                                               \
        }                                                               \
        Scalar* mdata() { return x; }                                   \
        const Scalar* mdata() const { return x; }                       \
        int mrows() const { return M; }                                 \
        int mcols() const { return N; }                                 \
        int mrow_stride() const { return RS; }                          \
        int mcol_stride() const { return CS; }                          \
        Scalar at(int i, int j) const { return x[i*RS+j*CS]; }          \
        Scalar& at(int i, int j) { return x[i*RS+j*CS]; }               \
        const row_t row(int i) const { return const_cast<self_t&>(*this).row(i); } \
        const col_t col(int i) const { return const_cast<self_t&>(*this).col(i); } \
        transpose_t transpose() { return transpose_t(*this); }          \
        const transpose_t transpose() const { return transpose_t(const_cast<self_t&>(*this)); }                


    //
    // Static rows, static cols
    //
    
    template <int M, int N, int RS, int CS, class Scalar>
    class RefMatrix : public FixedMatrix<M,N,RefMatrix<M,N,RS,CS,Scalar> >
    {
    private:
        Scalar *x;
    public:
        typedef RefMatrix<M,N,RS,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        row_t row(int i) { return row_t(x + i*RS); }
        col_t col(int i) { return col_t(x + i*CS); }
    };

    template <int M, int N, int CS, class Scalar>
    class RefMatrix<M,N,-1,CS,Scalar> : public FixedMatrix<M,N,RefMatrix<M,N,-1,CS,Scalar> >
    {
    private:
        Scalar *x;
        int RS;
    public:
        typedef RefMatrix<M,N,-1,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int RS_) : x(x_), RS(RS_) {}
        row_t row(int i) { return row_t(x + i*RS); }
        col_t col(int i) { return col_t(x + i*CS, RS); }
    };
    
    template <int M, int N, int RS, class Scalar>
    class RefMatrix<M,N,RS,-1,Scalar> : public FixedMatrix<M,N,RefMatrix<M,N,RS,-1,Scalar> >
    {
    private:
        Scalar *x;
        int CS;
    public:
        typedef RefMatrix<M,N,RS,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int CS_) : x(x_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, CS); }
        col_t col(int i) { return col_t(x + i*CS); }
    };

    template <int M, int N, class Scalar>
    class RefMatrix<M,N,-1,-1,Scalar> : public FixedMatrix<M,N,RefMatrix<M,N,-1,-1,Scalar> >
    {
    private:
        Scalar *x;
        int RS, CS;
    public:
        typedef RefMatrix<M,N,-1,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int RS_, int CS_) : x(x_), RS(RS_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, CS); }
        col_t col(int i) { return col_t(x + i*CS, RS); }
    };

    //
    // Dynamic rows, static cols
    //

    template <int N, int RS, int CS, class Scalar>
    class RefMatrix<-1,N,RS,CS,Scalar> : public FixedWidthMatrix<N,AbstractMatrix<RefMatrix<-1,N,RS,CS,Scalar> > >
    {
    private:
        Scalar *x;
        int M;
    public:
        typedef RefMatrix<-1,N,RS,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_) : x(x_), M(M_) {}
        row_t row(int i) { return row_t(x + i*RS); }
        col_t col(int i) { return col_t(x + i*CS, M); }
    };

    template <int N, int CS, class Scalar>
    class RefMatrix<-1,N,-1,CS,Scalar> : public FixedWidthMatrix<N,AbstractMatrix<RefMatrix<-1,N,-1,CS,Scalar> > >
    {
    private:
        Scalar *x;
        int M, RS;
    public:
        typedef RefMatrix<-1,N,-1,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int RS_) : x(x_), M(M_), RS(RS_) {}
        row_t row(int i) { return row_t(x + i*RS); }
        col_t col(int i) { return col_t(x + i*CS, M, RS); }
    };
    
    template <int N, int RS, class Scalar>
    class RefMatrix<-1,N,RS,-1,Scalar> : public FixedWidthMatrix<N,AbstractMatrix<RefMatrix<-1,N,RS,-1,Scalar> > >
    {
    private:
        Scalar *x;
        int M, CS;
    public:
        typedef RefMatrix<-1,N,RS,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int CS_) : x(x_), M(M_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, CS); }
        col_t col(int i) { return col_t(x + i*CS, M); }
    };

    template <int N, class Scalar>
    class RefMatrix<-1,N,-1,-1,Scalar> : public FixedWidthMatrix<N,AbstractMatrix<RefMatrix<-1,N,-1,-1,Scalar> > >
    {
    private:
        Scalar *x;
        int M, RS, CS;
    public:
        typedef RefMatrix<-1,N,-1,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int RS_, int CS_) : x(x_), M(M_), RS(RS_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, CS); }
        col_t col(int i) { return col_t(x + i*CS, M, RS); }
    };

    //
    // Static rows, dynamic cols
    //

    template <int M, int RS, int CS, class Scalar>
    class RefMatrix<M,-1,RS,CS,Scalar> : public FixedHeightMatrix<M,AbstractMatrix<RefMatrix<M,-1,RS,CS,Scalar> > >
    {
    private:
        Scalar *x;
        int N;
    public:
        typedef RefMatrix<M,-1,RS,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int N_) : x(x_), N(N_) {}
        row_t row(int i) { return row_t(x + i*RS, N); }
        col_t col(int i) { return col_t(x + i*CS); }
    };

    template <int M, int CS, class Scalar>
    class RefMatrix<M,-1,-1,CS,Scalar> : public FixedHeightMatrix<M,AbstractMatrix<RefMatrix<M,-1,-1,CS,Scalar> > >
    {
    private:
        Scalar *x;
        int N, RS;
    public:
        typedef RefMatrix<M,-1,-1,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int N_, int RS_) : x(x_), N(N_), RS(RS_) {}
        row_t row(int i) { return row_t(x + i*RS, N); }
        col_t col(int i) { return col_t(x + i*CS, RS); }
    };
    
    template <int M, int RS, class Scalar>
    class RefMatrix<M,-1,RS,-1,Scalar> : public FixedHeightMatrix<M,AbstractMatrix<RefMatrix<M,-1,RS,-1,Scalar> > >
    {
    private:
        Scalar *x;
        int N, CS;
    public:
        typedef RefMatrix<M,-1,RS,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int N_, int CS_) : x(x_), N(N_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, N, CS); }
        col_t col(int i) { return col_t(x + i*CS); }
    };

    template <int M, class Scalar>
    class RefMatrix<M,-1,-1,-1,Scalar> : public FixedHeightMatrix<M,AbstractMatrix<RefMatrix<M,-1,-1,-1,Scalar> > >
    {
    private:
        Scalar *x;
        int N, RS, CS;
    public:
        typedef RefMatrix<M,-1,-1,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int N_, int RS_, int CS_) : x(x_), N(N_), RS(RS_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, N, CS); }
        col_t col(int i) { return col_t(x + i*CS, RS); }
    };
    

    //
    // Dynamic rows, dynamic cols
    //

    template <int RS, int CS, class Scalar>
    class RefMatrix<-1,-1,RS,CS,Scalar> : public DynamicMatrix<RefMatrix<-1,-1,RS,CS,Scalar> >
    {
    private:
        Scalar *x;
        int M, N;
    public:
        typedef RefMatrix<-1,-1,RS,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int N_) : x(x_), M(M_), N(N_) {}
        row_t row(int i) { return row_t(x + i*RS, N); }
        col_t col(int i) { return col_t(x + i*CS, M); }
    };

    template <int CS, class Scalar>
    class RefMatrix<-1,-1,-1,CS,Scalar> : public DynamicMatrix<RefMatrix<-1,-1,-1,CS,Scalar> >
    {
    private:
        Scalar *x;
        int M, N, RS;
    public:
        typedef RefMatrix<-1,-1,-1,CS,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int N_, int RS_) : x(x_), M(M_), N(N_), RS(RS_) {}
        row_t row(int i) { return row_t(x + i*RS, N); }
        col_t col(int i) { return col_t(x + i*CS, M, RS); }
    };
    
    template <int RS, class Scalar>
    class RefMatrix<-1,-1,RS,-1,Scalar> : public DynamicMatrix<RefMatrix<-1,-1,RS,-1,Scalar> >
    {
    private:
        Scalar *x;
        int M, N, CS;
    public:
        typedef RefMatrix<-1,-1,RS,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int N_, int CS_) : x(x_), M(M_), N(N_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, N, CS); }
        col_t col(int i) { return col_t(x + i*CS, M); }
    };

    template <class Scalar>
    class RefMatrix<-1,-1,-1,-1,Scalar> : public DynamicMatrix<RefMatrix<-1,-1,-1,-1,Scalar> >
    {
    private:
        Scalar *x;
        int M, N, RS, CS;
    public:
        typedef RefMatrix<-1,-1,-1,-1,Scalar> self_t;
        LATL_REFM_COMMON();
        RefMatrix(Scalar *x_, int M_, int N_, int RS_, int CS_) : x(x_), M(M_), N(N_), RS(RS_), CS(CS_) {}
        row_t row(int i) { return row_t(x + i*RS, N, CS); }
        col_t col(int i) { return col_t(x + i*CS, M, RS); }
    };

#undef LATL_REFM_COMMON

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

        template <class T>
        Matrix(const AbstractMatrix<T>& other) {
            assign(other);
        }

        Matrix(const Matrix& other) {
            *this = other;
        }

        template <class E>
        Matrix(const MatrixExpr<E>& e) {
            e.instance()(*this);
        }

        Matrix& operator=(const Matrix& other) {
            assign(other);
            return *this;
        }

        template <class T>
        Matrix& operator=(const AbstractMatrix<T>& other) {
            assign(other);
            return *this;
        }
        
        
        Scalar* mdata() { return x; }
        const Scalar* mdata() const { return x; }

        int mrow_stride() const { return N; }
        int mcol_stride() const { return 1; }
        
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
    
    template <int N, class Scalar>
    class Matrix<-1,N,Scalar> : public FixedWidthMatrix<N,AbstractMatrix<Matrix<-1,N,Scalar> > >
    {
    private:
        Scalar *x;
        int M;

        void init(int m) {
            M = m;
            x = new Scalar[M*N];
        }
    public:
        typedef Matrix<-1,N,Scalar> self_t;
        typedef matrix_traits<self_t> traits;
        typedef typename traits::row_t row_t;
        typedef typename traits::col_t col_t;
        typedef typename traits::transpose_t transpose_t;

        Matrix() { M = 0; x = 0; }
        Matrix(int m) { init(m); }        
        Matrix(int m, Scalar s) {
            init(m);
            fill(*this, s);
        }

        void resize(int m) {
            Scalar *y = new Scalar[m*N];
            for (int i=0; i<(m<M?m:M); ++i)
                for (int j=0; j<N; ++j)
                    y[i*N+j] = x[i*N+j];            
            delete[] x;
            M = m;
            x = y;
        }

        template <class T>
        Matrix(const AbstractMatrix<T>& other) {
            init(other.rows());
            assign(other);
        }
        
        Matrix(const Matrix& other) {
            init(other.rows());
            assign(other);
        }

        template <class E>
        Matrix(const MatrixExpr<E>& e) {
            init(e.instance().rows());
            e.instance()(*this);
        }        

        Matrix& operator=(const Matrix& other) { assign(other); return *this; }
        
        template <class T>
        Matrix& operator=(const AbstractMatrix<T>& other) {
            assign(other);
            return *this;
        }
        
        int mrows() const { return M; }

        Scalar* mdata() { return x; }
        const Scalar* mdata() const { return x; }

        int mrow_stride() const { return N; }
        int mcol_stride() const { return 1; }
        
        Scalar at(int i, int j) const { return x[i*N + j]; }
        Scalar& at(int i, int j) { return x[i*N + j]; }

        row_t row(int i) { return row_t(x+i*N); }
        const row_t row(int i) const { return const_cast<self_t&>(*this).row(i); }
        
        col_t col(int i) { return col_t(x+i, M); }
        const col_t col(int i) const { return const_cast<self_t&>(*this).col(i); }

        transpose_t transpose() { return transpose_t(*this); }        
        const transpose_t transpose() const { return transpose_t(const_cast<self_t&>(*this)); }
    };

    template <int M, class Scalar>
    class Matrix<M,-1,Scalar> : public FixedHeightMatrix<M,AbstractMatrix<Matrix<M,-1,Scalar> > >
    {
    private:
        Scalar *x;
        int N;

        void init(int n) {
            N = n;
            x = new Scalar[M*N];
        }
    public:
        typedef Matrix<M,-1,Scalar> self_t;
        typedef matrix_traits<self_t> traits;
        typedef typename traits::row_t row_t;
        typedef typename traits::col_t col_t;
        typedef typename traits::transpose_t transpose_t;

        Matrix() { N = 0; x = 0; }
        Matrix(int n) { init(n); }        
        Matrix(int n, Scalar s) {
            init(n);
            fill(*this, s);
        }

        void resize(int n) {
            Scalar *y = new Scalar[M*n];
            for (int i=0; i<M; ++i)
                for (int j=0; j<(n<N?n:N); ++j)
                    y[i*n+j] = x[i*N+j];            
            delete[] x;
            N = n;
            x = y;
        }

        template <class T>
        Matrix(const AbstractMatrix<T>& other) {
            init(other.cols());
            assign(other);
        }
        
        Matrix(const Matrix& other) {
            init(other.cols());
            *this = other;
        }

        template <class E>
        Matrix(const MatrixExpr<E>& e) {
            init(e.instance().cols());
            e.instance()(*this);
        }        
        
        Matrix& operator=(const Matrix& other) { assign(other); return *this; }
        template <class T>
        Matrix& operator=(const AbstractMatrix<T>& other) {
            assign(other);
            return *this;
        }

        int mcols() const { return N; }

        Scalar* mdata() { return x; }
        const Scalar* mdata() const { return x; }

        int mrow_stride() const { return N; }
        int mcol_stride() const { return 1; }
        
        Scalar at(int i, int j) const { return x[i*N + j]; }
        Scalar& at(int i, int j) { return x[i*N + j]; }

        row_t row(int i) { return row_t(x+i*N, N); }
        const row_t row(int i) const { return const_cast<self_t&>(*this).row(i); }
        
        col_t col(int i) { return col_t(x+i, N); }
        const col_t col(int i) const { return const_cast<self_t&>(*this).col(i); }

        transpose_t transpose() { return transpose_t(*this); }        
        const transpose_t transpose() const { return transpose_t(const_cast<self_t&>(*this)); }
    };
    
    template <class Scalar>
    class Matrix<-1,-1,Scalar> : public DynamicMatrix<Matrix<-1,-1,Scalar> >
    {
    private:
        Scalar *x;
        int M,N;

        void init(int m, int n) {
            M = m;
            N = n;
            x = new Scalar[M*N];
        }
    public:
        typedef Matrix<-1,-1,Scalar> self_t;
        typedef matrix_traits<self_t> traits;
        typedef typename traits::row_t row_t;
        typedef typename traits::col_t col_t;
        typedef typename traits::transpose_t transpose_t;

        Matrix() { M = N = 0; x = 0; }
        Matrix(int m, int n) { init(m, n); }        
        Matrix(int m, int n, Scalar s) {
            init(m,n);
            fill(*this, s);
        }

        void resize(int m, int n) {
            Scalar *y = new Scalar[m*n];
            for (int i=0; i<(m<M?m:M); ++i)
                for (int j=0; j<(n<N?n:N); ++j)
                    y[i*n+j] = x[i*N+j];
            delete[] x;
            M = m;
            N = n;
            x = y;
        }

        template <class T>
        Matrix(const AbstractMatrix<T>& other) {
            init(other.rows(), other.cols());
            assign(other);
        }
        Matrix(const Matrix& other) {
            init(other.rows(), other.cols());
            assign(other);
        }

        template <class E>
        Matrix(const MatrixExpr<E>& e) {
            init(e.instance().rows(), e.instance().cols());
            e.instance()(*this);
        }        
        
        Matrix& operator=(const Matrix& other) { assign(other); return *this; }
        template <class T>
        Matrix& operator=(const AbstractMatrix<T>& other) {
            assign(other);
            return *this;
        }

        template <class E>
        Matrix& operator=(const MatrixExpr<E>& e) {
            e.instance()(*this);
            return *this;
        }        
        

        int mrows() const { return M; }
        int mcols() const { return N; }

        Scalar* mdata() { return x; }
        const Scalar* mdata() const { return x; }

        int mrow_stride() const { return N; }
        int mcol_stride() const { return 1; }
        
        Scalar at(int i, int j) const { return x[i*N + j]; }
        Scalar& at(int i, int j) { return x[i*N + j]; }

        row_t row(int i) { return row_t(x+i*N, N); }
        const row_t row(int i) const { return const_cast<self_t&>(*this).row(i); }
        
        col_t col(int i) { return col_t(x+i, M, N); }
        const col_t col(int i) const { return const_cast<self_t&>(*this).col(i); }

        transpose_t transpose() { return transpose_t(*this); }        
        const transpose_t transpose() const { return transpose_t(const_cast<self_t&>(*this)); }
    };

    //
    // Bounds checking for slices
    //

    template <int R0, int C0, int M, int N, class Mat>
    void check_slice_bounds(const AbstractMatrix<Mat>& m)
    {
        check_slice_bounds<R0,M>(m.T()[0]);
        check_slice_bounds<C0,N>(m[0]);
    }

    template <int M, int N, class Mat>
    void check_slice_bounds(const AbstractMatrix<Mat>& m, int R0, int C0)
    {
        check_slice_bounds<M>(m.T()[0], R0);
        check_slice_bounds<N>(m[0], C0);
    }                  

    //
    // Slices
    //

    template <int M, int N, int RS, int CS>
    struct MSliceCreator {
        template <class Mat> static
        RefMatrix<M,N,RS,CS, typename matrix_traits<Mat>::scalar_t>
        eval(AbstractMatrix<Mat>& m, int R0, int C0, int, int, int, int) {
            typedef typename matrix_traits<Mat>::scalar_t scalar_t;
            return RefMatrix<M,N,RS,CS,scalar_t>(m.data_at(R0,C0));
        }
    };

    template <int M, int N, int RS>
    struct MSliceCreator<M,N,RS,-1>
    {
        template <class Mat> static
        RefMatrix<M,N,RS,-1, typename matrix_traits<Mat>::scalar_t>
        eval(AbstractMatrix<Mat>& m, int R0, int C0, int, int, int, int CS) {
            typedef typename matrix_traits<Mat>::scalar_t scalar_t;
            return RefMatrix<M,N,RS,-1,scalar_t>(m.data_at(R0,C0), CS);
        }
    };    

    template <int M, int N, int CS>
    struct MSliceCreator<M,N,-1,CS>
    {
        template <class Mat> static
        RefMatrix<M,N,-1,CS, typename matrix_traits<Mat>::scalar_t>
        eval(AbstractMatrix<Mat>& m, int R0, int C0, int, int, int RS, int) {
            typedef typename matrix_traits<Mat>::scalar_t scalar_t;
            return RefMatrix<M,N,-1,CS,scalar_t>(m.data_at(R0,C0), RS);
        }
    };        
    
    template <int M, int N>
    struct MSliceCreator<M,N,-1,-1>
    {
        template <class Mat> static
        RefMatrix<M,N,-1,-1, typename matrix_traits<Mat>::scalar_t>
        eval(AbstractMatrix<Mat>& m, int R0, int C0, int, int, int RS, int CS) {
            typedef typename matrix_traits<Mat>::scalar_t scalar_t;
            return RefMatrix<M,N,-1,-1,scalar_t>(m.data_at(R0,C0), RS, CS);
        }
    };    

    template <int R0, int C0, int M, int N, class Mat>
    RefMatrix<M,N,
              matrix_traits<Mat>::static_row_stride,
              matrix_traits<Mat>::static_col_stride,
              typename matrix_traits<Mat>::scalar_t>
    slice(AbstractMatrix<Mat>& m)
    {
        check_slice_bounds<R0,C0,M,N>(m);
        return MSliceCreator<M,N,
            matrix_traits<Mat>::static_row_stride,
            matrix_traits<Mat>::static_col_stride>::eval(m, R0, C0,
                                                         M, N,
                                                         m.row_stride(), m.col_stride());
    }

    template <int R0, int C0, int M, int N, class Mat>
    const RefMatrix<M,N,
                    matrix_traits<Mat>::static_row_stride,
                    matrix_traits<Mat>::static_col_stride,
                    typename matrix_traits<Mat>::scalar_t>
    slice(const AbstractMatrix<Mat>& m)
    {
        return slice<R0,C0,M,N>(const_cast<AbstractMatrix<Mat>&>(m));
    }

    template <int M, int N, class Mat>
    RefMatrix<M,N,
              matrix_traits<Mat>::static_row_stride,
              matrix_traits<Mat>::static_col_stride,
              typename matrix_traits<Mat>::scalar_t>
    slice(AbstractMatrix<Mat>& m, int R0, int C0)
    {
        check_slice_bounds<M,N>(m, R0, C0);
        return MSliceCreator<M,N,
            matrix_traits<Mat>::static_row_stride,
            matrix_traits<Mat>::static_col_stride>::eval(m, R0, C0,
                                                         M, N,
                                                         m.row_stride(), m.col_stride());
    }

    template <int M, int N, class Mat>
    const RefMatrix<M,N,
                    matrix_traits<Mat>::static_row_stride,
                    matrix_traits<Mat>::static_col_stride,
                    typename matrix_traits<Mat>::scalar_t>
    slice(const AbstractMatrix<Mat>& m, int R0, int C0)
    {
        return slice<-1,-1,M,N>(const_cast<AbstractMatrix<Mat>&>(m), R0, C0);
    }

    

    //
    // as_row(), as_col()
    //
    
    template <int N, class V>
    RefMatrix<N,1,vector_traits<V>::static_stride,vector_traits<V>::static_stride*N,
              typename vector_traits<V>::scalar_t>
    as_col(FixedVector<N,V>& v) {
        typedef typename vector_traits<V>::scalar_t scalar_t;
        const int ss = vector_traits<V>::static_stride;
        return RefMatrix<N,1,ss,ss*N,scalar_t>(v.data());
    }

    template <int N, class V>
    RefMatrix<1,N,vector_traits<V>::static_stride*N,vector_traits<V>::static_stride,
              typename vector_traits<V>::scalar_t>
    as_row(FixedVector<N,V>& v) {
        typedef typename vector_traits<V>::scalar_t scalar_t;
        const int ss = vector_traits<V>::static_stride;
        return RefMatrix<1,N,ss*N,ss,scalar_t>(v.data());
    }

    template <class V>
    RefMatrix<-1,1,vector_traits<V>::static_stride,-1,
              typename vector_traits<V>::scalar_t>
    as_col(DynamicVector<V>& v) {
        typedef typename vector_traits<V>::scalar_t scalar_t;
        const int ss = vector_traits<V>::static_stride;
        return RefMatrix<-1,1,ss,-1,scalar_t>(v.data(), v.size(), v.size() * ss);
    }

    template <class V>
    RefMatrix<1,-1,-1,vector_traits<V>::static_stride,
              typename vector_traits<V>::scalar_t>
    as_row(DynamicVector<V>& v) {
        typedef typename vector_traits<V>::scalar_t scalar_t;
        const int ss = vector_traits<V>::static_stride;
        return RefMatrix<1,-1,-1,ss,scalar_t>(v.data(), v.size(), v.stride() * v.size());
    }

    //
    // Container
    //

    template <int M, int N, class Scalar>
    struct MatrixHolder {
        typedef Matrix<M,N,Scalar> type;
        type value;
        MatrixHolder() {}
        template <class Mat>
        MatrixHolder(const AbstractMatrix<Mat>& m) : value(m) {}        
        MatrixHolder(int, int) {}
        void resize(int, int) {}
        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };

    template <int N, class Scalar>
    struct MatrixHolder<-1,N,Scalar> {
        typedef Matrix<-1,N,Scalar> type;
        type value;
        MatrixHolder() {}
        
        template <class Mat>
        MatrixHolder(const AbstractMatrix<Mat>& m) : value(m) {}
        MatrixHolder(int M, int) : value(M) {}
        void resize(int M, int) { value.resize(M); }
        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };

    template <int M, class Scalar>
    struct MatrixHolder<M,-1,Scalar> {
        typedef Matrix<M,-1,Scalar> type;
        type value;
        MatrixHolder() {}
        
        template <class Mat>
        MatrixHolder(const AbstractMatrix<Mat>& m) : value(m) {}
        MatrixHolder(int, int N) : value(N) {}
        void resize(int, int N) { value.resize(N); }
        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };    

    template <class Scalar>
    struct MatrixHolder<-1,-1,Scalar> {
        typedef Matrix<-1,-1,Scalar> type;
        type value;
        MatrixHolder() {}
        
        template <class Mat>
        MatrixHolder(const AbstractMatrix<Mat>& m) : value(m) {}
        MatrixHolder(int M, int N) : value(M,N) {}
        void resize(int M, int N) { value.resize(M,N); }
        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };    
    

    
}

#endif
