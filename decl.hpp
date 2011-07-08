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

#ifndef LATL_DECL_HPP
#define LATL_DECL_HPP

namespace latl {

    //
    // Vectors
    //
    
    template <class V>
    struct vector_traits;    
    
    template <class V>
    struct AbstractVector;

    template <int N, class V>
    struct FixedVector;

    template <class V>
    struct DynamicVector;

    template <int N=-1, class Scalar=double>
    class Vector;

    template <int N, int Stride, class Scalar>
    class RefVector;

    template <class E>
    struct VectorExpr;

    //
    // Matrices
    //
    
    template <class Mat>
    struct matrix_traits;

    template <class Mat>
    struct AbstractMatrix;

    template <int N, class Mat>
    struct FixedWidthMatrix;

    template <int M, class Mat>
    struct FixedHeightMatrix;

    template <int M, int N, class Mat>
    struct FixedMatrix;

    template <class Mat>
    struct DynamicMatrix;

    template <int M=-1, int N=M, class Scalar=double>
    class Matrix;

    template <class Mat>
    class TransposeMatrix;
    
    template <int M, int N, int RowStride, int ColStride, class Scalar>
    class RefMatrix;
    
    template <class E>
    struct MatrixExpr;    
}

#endif
