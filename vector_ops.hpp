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
#define LATL_WIDER_VS(V1,V2) typename Wider<LATL_VS(V1),LATL_VS(V2)>::type
    
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

        typedef Vector<MaxInt<vector_traits<V1>::static_size, vector_traits<V2>::static_size>::value,
                       LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorSum<V1,V2>::result_t
    operator+(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorSum<V1,V2>(a,b);
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

        typedef Vector<MaxInt<vector_traits<V1>::static_size, vector_traits<V2>::static_size>::value,
                       LATL_WIDER_VS(V1,V2)> result_t;
    };
    
    template <class V1, class V2>
    typename VectorDiff<V1,V2>::result_t
    operator-(const AbstractVector<V1>& a, const AbstractVector<V2>& b) {
        return VectorDiff<V1,V2>(a,b);
    }

}

#endif
