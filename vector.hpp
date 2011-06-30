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

#ifndef LATL_VECTOR_HPP
#define LATL_VECTOR_HPP

#include <latl/debug.hpp>
#include <latl/decl.hpp>

namespace latl
{

    //
    // AbstractVector definition
    //
    
    template <class V>
    struct AbstractVector {
	V& instance() { return static_cast<V&>(*this); }
	const V& instance() const { return static_cast<const V&>(*this); }

	typedef typename vector_traits<V>::scalar_t scalar_t;

	enum {
            static_size = vector_traits<V>::static_size,
            static_stride = vector_traits<V>::static_stride,
        };
	
	scalar_t operator[](int i) const {
            LATL_CHECK_BOUNDS(i, 0, size());
            return instance().at(i);
        }
	scalar_t& operator[](int i) {
            LATL_CHECK_BOUNDS(i, 0, size());
            return instance().at(i);
        }

	int size() const { return instance().vsize(); }
	int stride() const { return instance().vstride(); }
        
        scalar_t* data() { return instance().vdata(); }
        const scalar_t* data() const { return instance.vdata(); }

        template <class W>
        void assign(const AbstractVector<W>& other)
        {
            assert(size() == other.size());
            unchecked_assign(other);
        }

        template <class W>
        void unchecked_assign(const AbstractVector<W>& other)
        {
            for (int i=0; i<size(); ++i)
                (*this)[i] = other[i];
        }
    };

    // 
    // Fill for constructors
    // 
    
    template <class V, class S>
    void fill(AbstractVector<V>& v, S s)
    {
        for (int i=0; i<v.size(); ++i)
            v[i] = s;
    }

    //
    // FixedVector definition
    //
    
    template <int N, class V>
    struct FixedVector : public AbstractVector<V>
    {
        int vsize() const { return AbstractVector<V>::static_size; }        
        
        template <class W>
        V& operator=(const FixedVector<N,W>& other) {
            unchecked_assign(other);
            return this->instance();
        }        

        template <class W>
        V& operator=(const DynamicVector<W>& other) {
            assign(other);
            return this->instance();
        }

        template <class E>
        V& operator=(const VectorExpr<E>& e) {
            e.instance()(*this);
            return this->instance();
        }
    };

    //
    // DynamicVector definition
    //
    
    template <class V>
    struct DynamicVector : public AbstractVector<V>
    {
        template <class W>
        V& operator=(const AbstractVector<W>& other) {
            assign(other);
            return this->instance();
        }
        template <class E>
        V& operator=(const VectorExpr<E>& e) {
            e.instance()(*this);
            return this->instance();
        }
    };

    //
    // Traits for Vector class
    //

    template <int N, class Scalar>
    struct vector_traits<Vector<N,Scalar> >
    {
        typedef Scalar scalar_t;
        enum {
            static_size = N,
            static_stride = 1,
        };
    };

    template <class Scalar>
    struct vector_traits<Vector<-1,Scalar> >
    {
        typedef Scalar scalar_t;
        enum {
            static_size = -1,
            static_stride = 1,
        };
    };
    
    //
    // Standard fixed-size Vector
    //
    
    template <int N, class Scalar>
    class Vector : public FixedVector<N,Vector<N,Scalar> >
    {
    private:
        Scalar x[N];
    public:
        Vector() {}
        Vector(const Scalar * s) {
            for (int i=0; i<N; ++i)
                x[i] = s[i];
        }
        
        Vector(Scalar s) {
            fill(*this, s);
        }

        template <class V>
        Vector(const FixedVector<N,V>& v) {
            *this = v;
        }

        template <class V>
        Vector(const DynamicVector<V>& v) {
            *this = v;
        }

        template <class E>
        Vector(const VectorExpr<E>& e) {
            e.instance()(*this);
        }

        int vstride() const { return 1; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }
        
        Scalar at(int i) const { return x[i]; }
        Scalar& at(int i) { return x[i]; }

        using FixedVector<N,Vector<N,Scalar> >::operator=;
    };

    //
    // Standard dynamic-sized Vector
    //
    
    template <class Scalar>
    class Vector<-1,Scalar> : public DynamicVector<Vector<-1,Scalar> >
    {
    private:
        Scalar *x;
        int N;
        void init(int n) {
            N = n;
            x = new Scalar[N];
        }
    public:
        Vector() : x(0), N(0) {}
        Vector(int n) {
            init(n);
        }
        Vector(int n, Scalar s) {
            init(n);
            fill(*this, s);
        }

        template <int NV, class V>
        Vector(const FixedVector<NV,V>& v) {
            init(NV);
            *this = v;
        }

        template <class V>
        Vector(const DynamicVector<V>& v) {
            init(v.size());
            *this = v;
        }        

        template <class E>
        Vector(const VectorExpr<E>& e) {
            init(e.size());
            e.instance()(*this);
        }
        
        ~Vector() {
            delete[] x;
        }
        
        using DynamicVector<Vector<-1,Scalar> >::operator=;

        int vsize() const { return N; }
        int vstride() const { return 1; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }
        
        void resize(int n) {
            Scalar *y = new Scalar[n];
            int m = n < N ? n : N;
            for (int i=0; i<m; ++i)
                y[i] = x[i];
            delete[] x;
            x = y;
            N = n;
        }
        
        Scalar at(int i) const { return x[i]; }        
        Scalar& at(int i) { return x[i]; }
    };


    //
    // Traits for RefVectors
    //

    template <int N, int Stride, class Scalar>
    struct vector_traits<RefVector<N,Stride,Scalar> >
    {
        typedef Scalar scalar_t;
        enum {
            static_size = N,
            static_stride = Stride,
        };
    };


    //
    // RefVector definitions
    //

    // Static size and stride
    
    template <int N, int Stride, class Scalar>
    class RefVector : public FixedVector<N,RefVector<N,Stride,Scalar> >
    {
    private:
        Scalar *x;
    public:
        RefVector(Scalar *s) : x(s) {}
        using FixedVector<N,RefVector<N,Stride,Scalar> >::operator=;

        RefVector& operator=(const RefVector& other) {
            unchecked_assign(other);
            return *this;
        }

        int vstride() const { return Stride; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }

        typedef Scalar scalar_t;
        scalar_t at(int i) const { return x[i*Stride]; }
        scalar_t& at(int i) { return x[i*Stride]; }
    };
    
    // Static stride
    
    template <int Stride, class Scalar>
    class RefVector<-1,Stride,Scalar> : public DynamicVector<RefVector<-1,Stride,Scalar> >
    {
    private:
        Scalar *x;
        int N;
    public:
        RefVector(Scalar *s, int n) : x(s), N(n) {}        
        using DynamicVector<RefVector<-1,Stride,Scalar> >::operator=;

        RefVector& operator=(const RefVector& other) {
            assign(other);
            return *this;
        }
        
        int vsize() const { return N; }
        int vstride() const { return Stride; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }
        
        typedef Scalar scalar_t;
        scalar_t at(int i) const { return x[i*Stride]; }
        scalar_t& at(int i) { return x[i*Stride]; }
    };

    // Static size
    
    template <int N, class Scalar>
    class RefVector<N,-1,Scalar> : public FixedVector<N,RefVector<N,-1,Scalar> >
    {
    private:
        Scalar *x;
        int Stride;
    public:
        RefVector(Scalar *s, int stride) : x(s), Stride(stride) {}
        using FixedVector<N,RefVector<N,-1,Scalar> >::operator=;

        RefVector& operator=(const RefVector& other) {
            unchecked_assign(other);
            return *this;
        }

        int vstride() const { return Stride; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }        
        
        typedef Scalar scalar_t;
        scalar_t at(int i) const { return x[i*Stride]; }
        scalar_t& at(int i) { return x[i*Stride]; }
    };

    // Dynamic size and stride
    
    template <class Scalar>
    class RefVector<-1,-1,Scalar> : public DynamicVector<RefVector<-1,-1,Scalar> >
    {
    private:
        Scalar *x;
        int N, Stride;
    public:
        RefVector(Scalar *s, int n, int stride) : x(s), N(n), Stride(stride) {}
        using DynamicVector<RefVector<-1,-1,Scalar> >::operator=;

        RefVector& operator=(const RefVector& other) {
            assign(other);
            return *this;
        }
        
        int vsize() const { return N; }
        int vstride() const { return Stride; }
        Scalar* vdata() { return x; }
        const Scalar* vdata() const { return x; }        
        
        typedef Scalar scalar_t;
        scalar_t at(int i) const { return x[i*Stride]; }
        scalar_t& at(int i) { return x[i*Stride]; }
    };
    
    //
    // Bounds checking for slices
    //
    
    template <int Start, int Size, int N, class V>
    void check_slice_bounds(const FixedVector<N,V>& v)
    {
        Assert<(Start >= 0 && Size >= 0 && Start + Size <= N)>();
    }

    template <int Start, int Size, class V>
    void check_slice_bounds(const DynamicVector<V>& v)
    {
        assert(Start >= 0 && Size >= 0 && Start + Size <= v.size());
    }

    //
    // slice() definitions
    //

    template <int Size, int Stride>
    struct SliceCreator {
        template <class V> static
        RefVector<Size,Stride, typename V::scalar_t>
        eval(AbstractVector<V>& v,
             int Start, int, int)
        {
            return RefVector<Size,Stride, typename V::scalar_t>(v.data() + Start * Stride);
        }
    };
    
    template <int Stride>
    struct SliceCreator<-1,Stride> {
        template <class V> static
        RefVector<-1,Stride, typename V::scalar_t>
        eval(AbstractVector<V>& v, int Start, int Size, int)
        {
            return RefVector<-1,Stride, typename V::scalar_t>(v.data() + Start * Stride,
                                                              Size);
        }
    };

    template <int Size>
    struct SliceCreator<Size,-1> {
        template <class V> static
        RefVector<Size,-1, typename V::scalar_t>
        eval(AbstractVector<V>& v, int Start, int, int Stride)
        {
            return RefVector<Size,-1, typename V::scalar_t>(v.data() + Start * Stride,
                                                            Stride);
        }
    };

    template <>
    struct SliceCreator<-1,-1> {
        template <class V> static
        RefVector<-1,-1, typename V::scalar_t>
        eval(AbstractVector<V>& v, int Start, int Size, int Stride)
        {
            return RefVector<-1,-1, typename V::scalar_t>(v.data() + Start * Stride,
                                                          Size,
                                                          Stride);
        }
    };
    
    template <int Start, int Size, class V>
    RefVector<Size,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(AbstractVector<V>& v)
    {
        check_slice_bounds<Start,Size>(v.instance());
        return SliceCreator<Size,
            vector_traits<V>::static_stride>::eval(v, Start, Size, v.stride());
    }

    template <int Start, int Size, class V>
    const RefVector<Size,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(const AbstractVector<V>& v)
    {
        check_slice_bounds<Start,Size>(v.instance());
        return SliceCreator<Size,
            vector_traits<V>::static_stride>::eval(const_cast<AbstractVector<V>&>(v),
                                                   Start, Size, v.stride());
    }

    template <class V>
    RefVector<-1,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(int Start, int Size, AbstractVector<V>& v)
    {
        assert(Start >= 0 && Size >= 0 && Start + Size <= v.size());
        return SliceCreator<-1,
            vector_traits<V>::static_stride>::eval(v, Start, Size, v.stride());
    }

    template <class V>
    const RefVector<-1,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(int Start, int Size, const AbstractVector<V>& v)
    {
        assert(Start >= 0 && Size >= 0 && Start + Size <= v.size());
        return SliceCreator<-1,
            vector_traits<V>::static_stride>::eval(const_cast<AbstractVector<V>&>(v),
                                                   Start, Size, v.stride());
    }
};

#endif
