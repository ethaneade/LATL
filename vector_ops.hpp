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

#ifndef LATL_VECTOR_OPS_HPP
#define LATL_VECTOR_OPS_HPP

#include <latl/debug.hpp>
#include <latl/vector.hpp>
#include <latl/expr.hpp>
#include <latl/scalar.hpp>

namespace latl
{
#define LATL_VS(V) typename vector_traits<V>::scalar_t
#define LATL_WIDER_VS(V1,V2) typename Wider< LATL_VS(V1) , LATL_VS(V2) >::type
    
    template <class V1, class V2>
    LATL_WIDER_VS(V1,V2) inner_product(const AbstractVector<V1>& a, const AbstractVector<V2>& b)
    {
        assert_same_size(a.instance(), b.instance());
        LATL_WIDER_VS(V1,V2) sum = 0;
        for (int i=0; i<a.size(); ++i)
            sum += a[i] * b[i];
        return sum;
    }

    template <class V1, class V2>
    LATL_WIDER_VS(V1,V2) operator*(const AbstractVector<V1>& a, const AbstractVector<V2>& b)
    {
        return inner_product(a,b);
    }

    // Scalar * Vector
    
    template <class S, class V>
    struct ScalarVectorProduct : public VectorExpr<ScalarVectorProduct<S,V> > {
        const AbstractVector<V>& v;
        S s;

        ScalarVectorProduct(const AbstractVector<V>& v_, S s_) : v(v_), s(s_) {}
        int size() const { return v.size(); }
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, v);
            for (int i=0; i<v.size(); ++i)
                w[i] = v[i] * s;
        }

        typedef Vector<vector_traits<V>::static_size,
                       typename Wider<S,LATL_VS(V)>::type> result_t;
    };
   
    template <class S, class V>
    typename ScalarVectorProduct<typename ScalarType<S>::type, V>::result_t
    operator*(S s, const AbstractVector<V>& v) {
        return ScalarVectorProduct<S,V>(v,s);
    }

    template <class S, class V>
    typename ScalarVectorProduct<typename ScalarType<S>::type, V>::result_t
    operator*(const AbstractVector<V>& v, S s) {
        return ScalarVectorProduct<S,V>(v,s);
    }

    template <class V> template <class S>
    V& AbstractVector<V>::operator*=(S s) {
        for (int i=0; i<size(); ++i)
            (*this)[i] *= s;
        return instance();
    }
    

    //
    // Vector / Scalar 
    //

    
    // For invertible types, multiply by the reciprocal
    template <class S, class V>    
    typename ScalarVectorProduct<typename ScalarType<S>::recip_type,V>::result_t
    operator/(const AbstractVector<V>& v, S s) {        
        return ScalarVectorProduct<S,V>(v,1/s);
    }    

    template <class V, class S>
    struct VectorQuotient : public VectorExpr<VectorQuotient<V,S> > {
        const AbstractVector<V>& v;
        S s;

        VectorQuotient(const AbstractVector<V>& v_, S s_) : v(v_), s(s_) {}
        int size() const { return v.size(); }
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, v);
            for (int i=0; i<v.size(); ++i)
                w[i] = v[i]/s;
        }

        typedef Vector<vector_traits<V>::static_size,
                       typename Wider<S,LATL_VS(V)>::type> result_t;
    };

    // For non-invertible types, divide each element
    template <class S, class V>
    typename VectorQuotient<V,typename ScalarType<S>::nonrecip_type>::result_t
    operator/(const AbstractVector<V>& v, S s) {
        return VectorQuotient<V,S>(v,s);
    }    

    template <class W, class V>
    void divide(AbstractVector<V>& v, typename ScalarType<W>::recip_type s) {
        W r = 1/s;
        for (int i=0; i<v.size(); ++i)
            v[i] *= r;
    }

    template <class W, class V>
    void divide(AbstractVector<V>& v, typename ScalarType<W>::nonrecip_type s) {
        for (int i=0; i<v.size(); ++i)
            v[i] /= s;
    }
    
    template <class V> template <class S>
    V& AbstractVector<V>::operator/=(S s) {
        typedef typename Wider<LATL_VS(V),S>::type W;
        divide<W>(*this, s);
        return instance();
    }

    
    
    //
    // Vector negation
    //
     
    template <class V>
    struct VectorNegation : public VectorExpr<VectorNegation<V> > {
        const AbstractVector<V>& v;
        int size() const { return v.size(); }
        VectorNegation(const AbstractVector<V>& v_) : v(v_) {}
        
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, v);
            for (int i=0; i<v.size(); ++i)
                w[i] = -v[i];
        }
        typedef Vector<vector_traits<V>::static_size,LATL_VS(V)> result_t;
    };
    
    template <class V>
    typename VectorNegation<V>::result_t
    operator-(const AbstractVector<V>& v) {
        return VectorNegation<V>(v);
    }

    //
    // Vector sum
    //
     
    template <class V1, class V2>
    struct VectorSum : public VectorExpr<VectorSum<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        int size() const { return a.size(); }
        VectorSum(const AbstractVector<V1>& a_,
                  const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
            assert_same_size(a,b);
        }
        
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, a);
            for (int i=0; i<a.size(); ++i)
                w[i] = a[i] + b[i];
        }

        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = MaxInt<N1,N2>::value,
        };

        typedef Vector<N,LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorSum<V1,V2>::result_t
    operator+(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorSum<V1,V2>(a,b);
    }

    template <class V1> template <class V2>
    V1& AbstractVector<V1>::operator+=(const AbstractVector<V2>& v) {
        assert_same_size(*this, v);
        for (int i=0; i<size(); ++i)
            (*this)[i] += v[i];
        return instance();
    }
    
    //
    // Vector difference
    //
     
    template <class V1, class V2>
    struct VectorDiff : public VectorExpr<VectorDiff<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        int size() const { return a.size(); }
        VectorDiff(const AbstractVector<V1>& a_,
                   const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
            assert_same_size(a,b);
        }
        
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, a);
            for (int i=0; i<a.size(); ++i)
                w[i] = a[i] - b[i];
        }

        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = MaxInt<N1,N2>::value,
        };

        typedef Vector<N,LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorDiff<V1,V2>::result_t
    operator-(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorDiff<V1,V2>(a,b);
    }

    template <class V1> template <class V2>
    V1& AbstractVector<V1>::operator-=(const AbstractVector<V2>& v) {
        assert_same_size(*this, v);
        for (int i=0; i<size(); ++i)
            (*this)[i] -= v[i];
        return instance();
    }

    //
    // Cross Product
    //    
    
    template <class V1, class V2>
    Vector<3,LATL_WIDER_VS(V1,V2)>
    operator^(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        assert_size_is<3>(a);
        assert_size_is<3>(b);
        Vector<3,LATL_WIDER_VS(V1,V2)> w;
        w[0] = a[1]*b[2] - a[2]*b[1];
        w[1] = a[2]*b[0] - a[0]*b[2];
        w[2] = a[0]*b[1] - a[1]*b[0];
        return w;
    }

    //
    // Vector concatenation
    //

    template <class V1, class V2>
    struct VectorConcat : public VectorExpr<VectorConcat<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = (N1 == -1 || N2 == -1) ? -1 : (N1+N2)
        };
        int size() const { return a.size() + b.size(); }
        VectorConcat(const AbstractVector<V1>& a_,
                   const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
        }

        template <class W>
        void operator()(AbstractVector<W>& w) const {
            CheckEquality<vector_traits<W>::static_size, N>::eval(w.size(), size());

            for (int i=0; i<a.size(); ++i)
                w[i] = a[i];
            for (int i=0; i<b.size(); ++i)
                w[i + a.size()] = b[i];
        }

        typedef Vector<N, LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorConcat<V1,V2>::result_t
    operator,(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorConcat<V1,V2>(a,b);
    }


    //
    // Vector element-wise multiplication
    //
     
    template <class V1, class V2>
    struct VectorDiagMult : public VectorExpr<VectorDiagMult<V1,V2> > {
        const AbstractVector<V1>& a;
        const AbstractVector<V2>& b;
        int size() const { return a.size(); }
        VectorDiagMult(const AbstractVector<V1>& a_,
                       const AbstractVector<V2>& b_)
            : a(a_), b(b_)
        {
            assert_same_size(a,b);
        }
        
        template <class W>
        void operator()(AbstractVector<W>& w) const {
            assert_same_size(w, a);
            for (int i=0; i<a.size(); ++i)
                w[i] = a[i] * b[i];
        }

        enum {
            N1 = vector_traits<V1>::static_size,
            N2 = vector_traits<V2>::static_size,
            N = MaxInt<N1,N2>::value,
        };

        typedef Vector<N,LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorDiagMult<V1,V2>::result_t
    diagmult(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorDiagMult<V1,V2>(a,b);
    }
    
}

#endif
