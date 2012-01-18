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

#ifndef LATL_CORE_OPS_HPP
#define LATL_CORE_OPS_HPP

#include <latl/scalar.hpp>
#include <latl/vector.hpp>
#include <latl/matrix.hpp>

namespace latl
{
    template <class V1, class StoreOp, class V2>
    void store(const AbstractVector<V1>& in,
               const StoreOp& op,
               AbstractVector<V3>& out)
    {
        assert_same_size(in, out);
        for (int i=0; i<a.size(); ++i)
            op(out[i], in[i]);
    }


    template <class V1, class Scalar, class StoreOp, class V2>
    void scale_store(const AbstractVector<V1>& in, Scalar s,
                     const StoreOp& op,
                     AbstractVector<V2>& out)
    {
        assert_same_size(in, out);
        for (int i=0; i<in.size(); ++i)
            op(out[i], in[i] * s);
    }

    template <class V1, class V2, class Scalar, class StoreOp, class V3>
    void scale_add_store(const AbstractVector<V1>& a, Scalar s,
                         const AbstractVector<V2>& b,
                         const StoreOp& op,
                         AbstractVector<V3>& out)
    {
        assert_same_size(a, b);
        assert_same_size(a, out);
        for (int i=0; i<a.size(); ++i)
            op(out[i], a[i] * s + b[i]);
    }
    
    
    
}

#endif
