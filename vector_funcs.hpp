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

#ifndef LATL_VECTOR_FUNCS_HPP
#define LATL_VECTOR_FUNCS_HPP

#include <latl/vector.hpp>
#include <latl/scalar.hpp>

namespace latl
{
    template <class S> Vector<1,typename ScalarType<S>::type>
    makeVector(S s0) {
        S s[] = { s0 };
        return Vector<1,typename ScalarType<S>::type>(s);
    }

    template <class S> Vector<2,typename ScalarType<S>::type>
    makeVector(S s0, S s1) {
        S s[] = { s0, s1 };
        return Vector<2,typename ScalarType<S>::type>(s);
    }
    
    template <class V>
    void zero(AbstractVector<V>& v) {
        fill(v, 0);
    }
    
    template <class V>
    LATL_VS(V) norm_sq(const AbstractVector<V>& v) { return v*v; }

    template <class V>
    LATL_VS(V) max_norm(const AbstractVector<V>& v) {
        LATL_VS(V) mn = 0;
        for (int i=0; i<v.size(); ++i)
            mn = latl::max(mn, latl::abs(v[i]));
        return mn;
    }
    
    template <class V>
    typename Wider<float,LATL_VS(V)>::type
    norm(const AbstractVector<V>& v)
    {
        return latl::sqrt(norm_sq(v));
    }
    
    template <class V>
    Vector<vector_traits<V>::static_size, LATL_VS(V)>
    unit(const AbstractVector<V>& v)
    {
        return v / norm(v);
    }

    template <class V>
    void normalize(AbstractVector<V>& v) {
        v /= norm(v);
    }


    template <int N, class V>
    Vector<N-1,LATL_VS(V)> project(const FixedVector<N,V>& v)
    {
        return slice<0,N-1>(v) / v[N-1];
    }

    template <int N, class V>
    Vector<N+1,LATL_VS(V)> append(const FixedVector<N,V>& v,
                                  LATL_VS(V) s)
    {
        Vector<N+1,LATL_VS(V)> vs;
        slice<0,N>(vs) = v;
        vs[N] = s;
        return vs;
    }

    template <int N, class V>
    Vector<N+1,LATL_VS(V)> unproject(const FixedVector<N,V>& v)
    {
        return append(v, (LATL_VS(V))1);
    }

    template <class V>
    Vector<2,LATL_VS(V)> perp(const FixedVector<2,V>& v)
    {
        Vector<2,LATL_VS(V)> pv;
        pv[0] = -v[1];
        pv[1] = v[0];
        return pv;
    }

    
    
}

#endif
