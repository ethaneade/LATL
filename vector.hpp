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
#include <latl/util.hpp>
#include <latl/vector_traits.hpp>

namespace latl
{

    template <bool DynN1, bool DynN2>
    struct AssertSameSize {
        template <class A, class B>
        static void eval(const A& a, const B& b) {
            assert(a.size() == b.size());
        }
    };

    template <>
    struct AssertSameSize<false,false> {
        template <class A, class B>
        static void eval(const A& a, const B& b) {
            Assert<((int)vector_traits<A>::static_size ==
                    (int)vector_traits<B>::static_size)>();
        }
    };

    template <class A, class B>
    void assert_same_size(const A& a, const B& b) {
        AssertSameSize<vector_traits<A>::static_size == -1,
            vector_traits<B>::static_size == -1>::eval(a,b);
    }

    template <int N, class V>
    void assert_size_is(const AbstractVector<V>& v) {
        CheckEquality<N,vector_traits<V>::static_size>::eval(N,v.size());
    }

    
    //
    // AbstractVector definition
    //
    
    template <class V>
    class AbstractVector {
    public:
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
            assert_same_size(instance(), other.instance());
            unchecked_assign(other);
        }

        template <class W>
        V& operator=(const AbstractVector<W>& w) {
            assign(w);
            return instance();
        }

        template <class W>
        void swap(AbstractVector<W>& w) {
            assert_same_size(instance(), w.instance());
            for (int i=0; i<size(); ++i) {
                scalar_t tmp = w[i];
                w[i] = (*this)[i];
                (*this)[i] = tmp;
            } 
        }

        template <class W> V& operator+=(const AbstractVector<W>& w);
        template <class W> V& operator-=(const AbstractVector<W>& w);
        template <class S> V& operator*=(S s);
        template <class S> V& operator/=(S s);
        
    private:
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

        using AbstractVector<V>::operator=;
        
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
        template <class E>
        V& operator=(const VectorExpr<E>& e) {
            e.instance()(*this);
            return this->instance();
        }
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

        Vector(const Vector& v) {
            assign(v);
        }
        
        template <class V>
        Vector(const AbstractVector<V>& v) {
            assign(v);
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

        Vector& operator=(const Vector& v) {
            assign(v);
            return *this;
        }

        template <class V>
        Vector& operator=(const AbstractVector<V>& v) {
            assign(v);
            return *this;
        }
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
        
        Vector(const Vector& v) {
            init(v.size());
            assign(v);
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
            init(e.instance().size());
            e.instance()(*this);
        }
        
        ~Vector() {
            delete[] x;
        }
        
        using DynamicVector<Vector<-1,Scalar> >::operator=;

        Vector& operator=(const Vector& v) {
            assign(v);
            return *this;
        }
        

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

        RefVector& operator=(const RefVector& other) {
            assign(other);
            return *this;
        }

        template <class W>
        RefVector& operator=(const AbstractVector<W>& other) {
            assign(other);
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
        template <class W>
        RefVector& operator=(const AbstractVector<W>& other) {
            assign(other);
            return *this;
        }

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
        
        template <class W>        
        RefVector& operator=(const AbstractVector<W>& other) {
            assign(other);
            return *this;
        }

        RefVector& operator=(const RefVector& other) {
            assign(other);
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

        template <class W>
        RefVector& operator=(const AbstractVector<W>& other) {
            assign(other);
            return *this;
        }

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
    // Bounds checking for slices:
    // Check as much as possible at compile time
    //

    template <int Start, int Size, int N>
    struct CheckSliceBounds {
        static void eval(int, int, int) {
            Assert<(Start >= 0 && Size >= 0 && Start + Size <= N)>();
        }
    };

    template <int Size, int N>
    struct CheckSliceBounds<-1,Size,N> {
        static void eval(int Start, int, int) {
            Assert<(Size >= 0 && Size <= N)>();
            assert(Start + Size <= N);
        }
    };

    template <int Start, int N>
    struct CheckSliceBounds<Start,-1,N> {
        static void eval(int, int Size, int) {
            Assert<(Start >= 0 && Start <= N)>();
            assert(Start + Size <= N);
        }
    };

    template <int N>
    struct CheckSliceBounds<-1,-1,N> {
        static void eval(int Start, int Size, int) {
            assert(Start >= 0 && Size >= 0 && Start + Size <= N);
        }
    };

    template <int Start, int Size>
    struct CheckSliceBounds<Start,Size,-1> {
        static void eval(int, int, int N) {
            Assert<(Start >= 0 && Size >= 0)>();
            assert(Start + Size <= N);
        }
    };

    template <int Size>
    struct CheckSliceBounds<-1,Size,-1> {
        static void eval(int Start, int, int N) {
            Assert<(Size >= 0)>();
            assert(Start >= 0 && Start + Size <= N);
        }
    };

    template <int Start>
    struct CheckSliceBounds<Start,-1,-1> {
        static void eval(int, int Size, int N) {
            Assert<(Start >= 0)>();
            assert(Size >= 0 && Start + Size <= N);
        }
    };
    
    template <>
    struct CheckSliceBounds<-1,-1,-1> {
        static void eval(int Start, int Size, int N) {
            assert(Start >= 0 && Size >= 0 && Start + Size <= N);
        }
    };
    

    template <int Start, int Size, class V>
    void check_slice_bounds(const AbstractVector<V>& v)
    {
        CheckSliceBounds<Start,Size,vector_traits<V>::static_size>::eval(Start,Size,v.size());
    }

    template <int Size, class V>
    void check_slice_bounds(const AbstractVector<V>& v, int Start)
    {
        CheckSliceBounds<-1,Size,vector_traits<V>::static_size>::eval(Start,Size,v.size());
    }
    
    //
    // slice() definitions
    //

    template <int Size, int Stride>
    struct VSliceCreator {
        template <class V> static
        RefVector<Size,Stride, typename V::scalar_t>
        eval(AbstractVector<V>& v,
             int Start, int, int)
        {
            return RefVector<Size,Stride, typename V::scalar_t>(v.data() + Start * Stride);
        }
    };
    
    template <int Stride>
    struct VSliceCreator<-1,Stride> {
        template <class V> static
        RefVector<-1,Stride, typename V::scalar_t>
        eval(AbstractVector<V>& v, int Start, int Size, int)
        {
            return RefVector<-1,Stride, typename V::scalar_t>(v.data() + Start * Stride,
                                                              Size);
        }
    };

    template <int Size>
    struct VSliceCreator<Size,-1> {
        template <class V> static
        RefVector<Size,-1, typename V::scalar_t>
        eval(AbstractVector<V>& v, int Start, int, int Stride)
        {
            return RefVector<Size,-1, typename V::scalar_t>(v.data() + Start * Stride,
                                                            Stride);
        }
    };

    template <>
    struct VSliceCreator<-1,-1> {
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
        return VSliceCreator<Size,
            vector_traits<V>::static_stride>::eval(v, Start, Size, v.stride());
    }

    template <int Start, int Size, class V>
    const RefVector<Size,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(const AbstractVector<V>& v)
    {
        check_slice_bounds<Start,Size>(v.instance());
        return VSliceCreator<Size,
            vector_traits<V>::static_stride>::eval(const_cast<AbstractVector<V>&>(v),
                                                   Start, Size, v.stride());
    }

    template <class V>
    RefVector<-1,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(int Start, int Size, AbstractVector<V>& v)
    {
        assert(Start >= 0 && Size >= 0 && Start + Size <= v.size());
        return VSliceCreator<-1,
            vector_traits<V>::static_stride>::eval(v, Start, Size, v.stride());
    }

    template <class V>
    const RefVector<-1,vector_traits<V>::static_stride,typename V::scalar_t>
    slice(int Start, int Size, const AbstractVector<V>& v)
    {
        assert(Start >= 0 && Size >= 0 && Start + Size <= v.size());
        return VSliceCreator<-1,
            vector_traits<V>::static_stride>::eval(const_cast<AbstractVector<V>&>(v),
                                                   Start, Size, v.stride());
    }

    //
    // Container
    //

    template <int N, class Scalar>
    struct VectorHolder {
        typedef Vector<N,Scalar> type;
        type value;
        VectorHolder() {}
        template <class V>
        VectorHolder(const AbstractVector<V>& v) : value(v) {}        
        VectorHolder(int) {}
        void resize(int) {}

        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };
 
    template <class Scalar>
    struct VectorHolder<-1,Scalar> {
        typedef Vector<-1,Scalar> type;
        type value;
        VectorHolder() {}
        
        template <class V>
        VectorHolder(const AbstractVector<V>& v) : value(v) {}
        VectorHolder(int N) : value(N) {}
        void resize(int N) { value.resize(N); }

        type& operator()() { return value; }
        const type& operator()() const { return value; }
    };    
};

#endif
