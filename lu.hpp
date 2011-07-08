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

#ifndef LATL_LU_HPP
#define LATL_LU_HPP

#include <latl/latl.hpp>

namespace latl
{
    template <int N, class Scalar>
    class LU
    {
    public:
        MatrixHolder<N,N,Scalar> lu;
        VectorHolder<N,Scalar> inv_diag;
        VectorHolder<N,int> index;
        int r;
        int swaps;
        
    public:
        LU() {}
        
        template <class Mat>
        LU(const AbstractMatrix<Mat>& m) { compute(m); }
        

        template <class Mat>
        void compute(const AbstractMatrix<Mat>& m) {
            assert_square(m);
            int n = m.rows();
            lu.resize(n, n);
            lu() = m;
            
            inv_diag.resize(n);
            index.resize(n);
            int row = 0;
            swaps = 0;

            for (int i=0; i<n; ++i)
                index()[i] = i;

            for (int i=0; i<n; ++i) {
                int argmax = row;
                Scalar maxval = latl::abs(lu.value(row,i));
                for (int j=row+1; j<n; ++j) {
                    Scalar ji = lu.value(j,i);
                    if (ji > maxval) {
                        maxval = ji;
                        argmax = j;
                    }
                }

                if (argmax != row) {
                    lu()[row].swap(lu()[argmax].instance());

                    int tmp = index()[i];
                    index()[i] = index()[argmax];
                    index()[argmax] = tmp;
                    ++swaps;
                }
                
                if (maxval < Constants<Scalar>::epsilon()) {
                    inv_diag()[i] = 0;
                    continue;
                }
                
                Scalar inv_ii = 1 / lu.value(row,i);
                inv_diag()[i] = inv_ii;
                
                for (int j=row+1; j<n; ++j) {
                    Scalar factor = lu.value(j,i) * inv_ii;
                    for (int k=i+1; k<n; ++k)
                        lu.value(j,k) -= factor * lu.value(row,k); 
                    lu.value(j,i) = factor;
                }
                ++row;                                               
            }
            r = row;
        }

        template <class V>
        typename VectorHolder<N,Scalar>::type
        inverse_times(const AbstractVector<V>& v) const {
            assert_same_size(index(), v);
            const int n = index().size();
            VectorHolder<N,Scalar> y(n);
            for (int i=0; i<n; ++i) {
                Scalar yi = v[index()[i]];
                for (int j=0; j<i; ++j)
                    yi -= lu.value(i, j) * y()[j];                    
                y()[i] = yi;
            }
            for (int i=n-1; i>=0; --i) {
                Scalar yi = y()[i];
                for (int j=i+1; j<n; ++j)
                    yi -= lu.value(i,j) * y()[j];
                y()[i] = inv_diag()[i] * yi;
            }
            return y.value;
        }

        int rank() const { return r; }

        bool is_full_rank() const { return r == index().size(); }
        
        Scalar determinant() const {
            Scalar det = lu()(0,0);
            for (int i=1; i<index().size(); ++i)
                det *= lu()(i,i);
            return ((swaps%2) ? -det : det);
        }
    };
}

#endif
