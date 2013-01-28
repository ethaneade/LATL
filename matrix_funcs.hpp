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
#include <latl/lu.hpp>

namespace latl
{
    template <class Mat>
    void zero(AbstractMatrix<Mat>& m) {
        fill(m, 0);
    }

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
    void symmetrize_from_upper(AbstractMatrix<Mat>& m)
    {
        assert_square(m);
        for (int i=0; i<m.rows(); ++i) {
            for (int j=0; j<i; ++j) {
                m(i,j) = m(j,i);
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
    LATL_MS(Mat) determinant(const FixedMatrix<2,2,Mat>& m)
    {
        return m(0,0)*m(1,1) - m(0,1)*m(1,0);
    }

    template <class Mat>
    LATL_MS(Mat) determinant(const FixedMatrix<3,3,Mat>& m)
    {
        LATL_MS(Mat) a, b, c;
        a = m(0,1)*m(1,2) - m(0,2)*m(1,1);
        b = m(0,0)*m(1,2) - m(0,2)*m(1,0);
        c = m(0,0)*m(1,1) - m(0,1)*m(1,0);
        return m(2,0)*a - m(2,1)*b + m(2,2)*c;
    }
    
    template <class Mat>
    Matrix<2,2,LATL_MS(Mat)> inverse(const FixedMatrix<2,2,Mat>& m)
    {
        typedef LATL_MS(Mat) Scalar;
        Scalar det = determinant(m);
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

    template <class M1, class M2, class Op>
    void outer_product_upper(const AbstractMatrix<M1>& J, LATL_MS(M1) s, const Op& op, AbstractMatrix<M2>& m)
    {
        assert_same_size(J.T()[0], m[0]);
        assert_square(m);
        for (int i=0; i<J.cols(); ++i)
            outer_product_upper(J.T()[i], s, op, m);
    }
    
    template <class Mat>
    LATL_MS(Mat) max_norm(const AbstractMatrix<Mat>& m) {
        LATL_MS(Mat) mn = 0;
        for (int i=0; i<m.rows(); ++i)
            mn = latl::max(mn, max_norm(m[i]));
        return mn;
    }

    //
    // Matrix exponential using scaling and squaring,
    // with diagonal Pade approximant for m=7.
    //
    template <class Mat1, class Mat2>
    void exp(const AbstractMatrix<Mat1>& m, AbstractMatrix<Mat2>& expm)
    {
        assert_square(m);
        assert_same_shape(m, expm);
        const int N = matrix_traits<Mat1>::static_rows;
        const int n = m.rows();
        
        typedef LATL_MS(Mat1) S;
        S scale = max_norm(m);
        S factor = 1;
        int s = 0;
        while (factor * scale > S(0.5)) {
            ++s;
            factor *= S(0.5);
        }
        MatrixHolder<N,N,S> sm(m * factor);
        MatrixHolder<N,N,S> term(sm());
        
        static const S c[7] = {S(1.0/2), S(3.0/26), S(5.0/312), S(5.0/3432),
                               S(1.0/11440), S(1.0/308880), S(1.0/17297280) };
        
        MatrixHolder<N,N,S> p(n,n);
        p() = Identity();
        MatrixHolder<N,N,S> q(p());
        
        S qfactor = -1;
        for (int i=0; i<7; ++i)
        {
            p() += c[i] * term();
            q() += qfactor * c[i] * term();
            term() = sm() * term();
            qfactor = -qfactor;
        }        
        LU<N,S> lu(q());
        assert(lu.is_full_rank());
        
        lu.inverse_times(p(), expm);
        
        for (int i=0; i<s; ++i)
            expm = expm * expm;
    }

    //
    // Compute S such that m == S*S.
    // Works for fixed or dynamic square matrices.
    // Returns false on failure to converge.
    //
    // Matrices with real negative eigenvalues have no real square root.
    //
    
    template <class Mat1, class Mat2>
    bool sqrt(const AbstractMatrix<Mat1>& m, AbstractMatrix<Mat2>& sqrt_m)
    {        
        typedef LATL_MS(Mat1) Scalar;
        assert_square(m);
        assert_same_shape(m, sqrt_m);

        const int N = matrix_traits<Mat1>::static_rows;
        const int n = m.rows();

        // Denman-Beavers iteration with scaling, as described in:
        //    Nicholas J. Higham. "Stable Iterations for the Matrix Square Root".
        //    Numerical Algorithms 15, 1997.
        //
        // lambda_k = |(det(Y_k) * det(Z_k))| ^ {-1/(2N)}
        // Y_{k+1} = 1/2 * (lambda_k * Y_k + 1/lambda_k * Z_k^{-1})
        // Z_{k+1} = 1/2 * (lambda_k * Z_k + 1/lambda_k * Y_k^{-1})
        //        
        // with Y_0 = m, Z_0 = I.
        // sqrt_m stands in for Y in the code below, as Y --> m^{1/2}.
        //
        
        // Scaling exponent:
        const Scalar lambda_p = -Scalar(1) / (2*N);
        
        MatrixHolder<N,N,Scalar> Z(n,n), Yinv(n,n), Zinv(n,n);        
        LU<N,Scalar> lu;

        // Y0, Z0
        sqrt_m = m;
        Z() = Identity();
        
        for (int iter=0; iter<32; ++iter)
        {
            Scalar detZ;
            
            if (iter) {
                lu.compute(Z());
                detZ = lu.determinant();
                lu.get_inverse(Zinv());
            } else {
                detZ = 1;
                Zinv() = Z();
            }
            
            lu.compute(sqrt_m);
            Scalar detY = lu.determinant();
            lu.get_inverse(Yinv());

            Scalar lambda = latl::pow(latl::abs(detY * detZ), lambda_p);
            Scalar inv_lambda = 1/lambda;
            
            // Y_{k+1} = 1/2 * (lambda_k * Y_k + 1/lambda_k * Z_k^{-1})
            sqrt_m *= lambda;
            Zinv() *= inv_lambda;
            sqrt_m += Zinv();
            sqrt_m *= Scalar(0.5);
            
            // Z_{k+1} = 1/2 * (lambda_k * Z_k + 1/lambda_k * Y_k^{-1})
            Z() *= lambda;
            Yinv() *= inv_lambda;
            Z() += Yinv();
            Z() *= Scalar(0.5);

            if ((iter % 4) == 3) {
                // At convergence, Y*Z = I.
                matrix_multiply(sqrt_m, Z(), ops::Assign(), Zinv());
                for (int i=0; i<n; ++i)
                    Zinv()(i,i) -= 1;
                if (max_norm(Zinv()) <= Constants<Scalar>::epsilon()*2)
                    return true;
            }
        }
        return false;
    }
    

    template <class Mat1, class Mat2>
    bool log(const AbstractMatrix<Mat1>& m, AbstractMatrix<Mat2>& log_m)
    {        
        typedef LATL_MS(Mat1) Scalar;
        assert_square(m);
        assert_same_shape(m, log_m);

        const int N = matrix_traits<Mat1>::static_rows;
        const int n = m.rows();
        
        MatrixHolder<N,N,Scalar> I(n,n);
        I() = Identity();

        Scalar scale = 1;
        MatrixHolder<N,N,Scalar> S(m), SS(n,n);
        while (max_norm(I() - S()) > Scalar(0.5))
        {            
            if (!sqrt(S(), SS()))
                return false;
            
            if (!sqrt(SS(), S()))
                return false;
            
            scale *= 4;
        }
        
        static const Scalar a[7] = {-1, 3, Scalar(-535.0/156),
                                    Scalar(145.0/78), Scalar(-1377.0/2860), Scalar(223.0/4290),
                                    Scalar(-11.0/7280)};
        
        static const Scalar b[7] = {Scalar(-7.0/2), Scalar(63.0/13), Scalar(-175.0/52),
                                    Scalar(175.0/143), Scalar(-63.0/286), Scalar(7.0/429),
                                    Scalar(-1.0/3432)};

        S() *= -1;
        for (int i=0; i<n; ++i)
            S()(i,i) += 1;
        
        MatrixHolder<N,N,Scalar> X(S());
        MatrixHolder<N,N,Scalar> P(n,n);
        zero(P());
        MatrixHolder<N,N,Scalar> Q(I());
        
        for (int i=0; i<7; ++i) {
            P() += a[i] * X();
            Q() += b[i] * X();
            X() = S() * X();
        }

        LU<N,Scalar> lu(Q());
        if (!lu.is_full_rank())
            return false;

        lu.inverse_times(P(), log_m);
       
        log_m *= scale;
        return true;
    }
    
}

#endif

