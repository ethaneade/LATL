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

#ifndef LATL_EXPR_HPP
#define LATL_EXPR_HPP

#include <latl/matrix.hpp>
#include <cstdlib>

namespace latl {

    template <class E>
    struct VectorExpr {
        E& instance() { return static_cast<E&>(*this); }
        const E& instance() const { return static_cast<const E&>(*this); }        
    };

    template <class E>
    struct MatrixExpr {
        E& instance() { return static_cast<E&>(*this); }
        const E& instance() const { return static_cast<const E&>(*this); }
    };
    
    struct Zeros : public VectorExpr<Zeros>, public MatrixExpr<Zeros> {
        template <class V>
        void operator()(AbstractVector<V>& v) const {
            fill(v, 0);
        }
        template <class M>
        void operator()(AbstractMatrix<M>& m) const {
            fill(m, 0);
        }
    };

    
    struct Random : public VectorExpr<Random>, public MatrixExpr<Random> {
        double factor;
        Random(double s = 1) : factor(s/RAND_MAX) {}
        template <class V>
        void operator()(AbstractVector<V>& v) const {
            for (int i=0; i<v.size(); ++i)
                v[i] = rand() * factor;
        }
        template <class M>
        void operator()(AbstractMatrix<M>& m) const {
            for (int i=0; i<m.rows(); ++i)
                for (int j=0; j<m.cols(); ++j)
                    m(i,j) = rand() * factor;
        }
    };

    struct CenteredRandom : public VectorExpr<CenteredRandom>, public MatrixExpr<CenteredRandom> {
        double a, b;
        CenteredRandom(double s = 1) : a(s / RAND_MAX), b(s*0.5) {}
        template <class V>
        void operator()(AbstractVector<V>& v) const {
            for (int i=0; i<v.size(); ++i)
                v[i] = rand() * a - b;
        }
        template <class M>
        void operator()(AbstractMatrix<M>& m) const {
            for (int i=0; i<m.rows(); ++i)
                for (int j=0; j<m.cols(); ++j)
                    m(i,j) = rand() * a - b;
        }
    };
    
    struct Identity : public MatrixExpr<Identity> {
        template <class M>
        void operator()(AbstractMatrix<M>& m) const {
            assert_square(m);
            fill(m, 0);
            for (int i=0; i<m.rows(); ++i)
                m(i,i) = 1;
        }
    };
    
    template <int N, class S, class Prev>
    struct VectorCreator : public VectorExpr<VectorCreator<N,S,Prev> >
    {
        S s;
        const Prev& prev;

        VectorCreator(S s_, const Prev& prev_) : s(s_), prev(prev_) {}
        int size() const { return N; }
        template <class V>
        void operator()(AbstractVector<V>& v) const {
            assert_size_is<N>(v);
            dump(v);
        }

        template <class V>
        void dump(AbstractVector<V>& v) const {
            v[N-1] = s;
            prev.dump(v);
        }
    };

    template <>
    struct VectorCreator<0,int,int> : public VectorExpr<VectorCreator<0,int,int> > {
        int size() const { return 0; }
        
        template <class V>
        void operator()(AbstractVector<V>& v) const {
            assert_size_is<0>(v);
        }

        template <class V>
        void dump(AbstractVector<V>& v) const {}
    };

    template <int N, class T, class S, class Prev>
    VectorCreator<N+1,T,VectorCreator<N,S,Prev> > operator,(const VectorCreator<N,S,Prev>& vc, T t)
    {
        return VectorCreator<N+1,T,VectorCreator<N,S,Prev> >(t, vc);
    }


    inline
    VectorCreator<0,int,int> vec() { return VectorCreator<0,int,int>(); }
    
}

#endif
