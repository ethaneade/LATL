#pragma once

#include <latl/ldlt.hpp>
#include <vector>

namespace latl {
    
    template <int N, typename S>
    struct BlockMatrix
    {
        Matrix<N,N,S> *m;
        int num_rows, num_cols;

        BlockMatrix(Matrix<N,N,S> *m_, int nr, int nc)
            : m(m_), num_rows(nr), num_cols(nc) {}

        int rows() const { return num_rows; }
        int cols() const { return num_cols; }

        template <class Mat>
        void get(int i, int j, FixedMatrix<N,N,Mat> &mij) const
        {
            mij = m[i*num_cols + j];
        }

        template <class Mat>
        void set(int i, int j, const FixedMatrix<N,N,Mat> &mij)
        {
            m[i*num_cols + j] = mij;
        }    
    };

    template <int N, typename S>
    class BlockLDLT
    {
    private:   
        std::vector<LDLT<N,S> > diag;
        std::vector<Matrix<N,N,S> > L;
        bool full_rank;
        bool pos_def;

    public:
        BlockLDLT() {
            full_rank = false;
            pos_def = false;
        }
        template <class BlockMat>
        void reconstitute(BlockMat &m) const
        {
            const int n = diag.size();
            for (int i=0; i<n; ++i) {
                for (int j=i; j<n; ++j) {
                    Matrix<N,N,S> mij = L[i*n+j];
                    for (int k=0; k<i; ++k) {
                        mij += L[i*n+k] * L[k*n + j];
                    }
                    m.set(i,j, mij);
                    m.set(j,i, mij.T());
                }
            }
        }

        bool is_full_rank() const { return full_rank; }
    
        bool is_positive_definite() const { return pos_def; }
    
        template <class BlockMat>
        bool compute(const BlockMat& m) {
            const int n = m.rows();
            assert(n == m.cols());

            full_rank = false;
            pos_def = true;
            L.resize(n*n);        
            diag.resize(n);
        
            for (int i=0; i<n; ++i) {
                Matrix<N,N,S> d;
                m.get(i,i, d);
                for (int k=0; k<i; ++k)
                    matrix_multiply(L[i*n+k], L[k*n+i], ops::Subtract(), d);

                if (!diag[i].compute(d)) {
                    pos_def = false;
                    return false;
                }
                pos_def = pos_def && diag[i].is_positive_definite();
            
                L[i*n+i] = d;
                
                for (int j=i+1; j<n; ++j) {
                    Matrix<N,N,S> y;
                    m.get(i,j,y);
                    for (int k=0; k<i; ++k)
                        matrix_multiply(L[i*n+k], L[k*n+j], ops::Subtract(), y);
                    diag[i].inverse_times(y, L[j*n+i].T().instance());
                    L[i*n+j] = y;
                }
            }
            full_rank = true;
            return true;
        }

        void inverse_times(const std::vector<Vector<N,S> > &x,
                           std::vector<Vector<N,S> > &y) const
        {
            assert(x.size() == diag.size());
            const int n = diag.size();
            y.resize(n);
        
            for (int i=0; i<n; ++i) {
                Vector<N,S> yi = x[i];
                for (int j=0; j<i; ++j)
                    matrix_vector_multiply(L[i*n+j], y[j], ops::Subtract(), yi);
                y[i] = yi;
            }
            for (int i=n-1; i>=0; --i) {
                Vector<N,S> yi = diag[i].inverse_times(y[i]);
                for (int j=i+1; j<n; ++j)
                    matrix_vector_multiply(L[j*n+i].T(), y[j], ops::Subtract(), yi);
                y[i] = yi;
            }
        }    
    };
}


#if LATL_BLOCK_LDLT_TEST
#include "io.hpp"
using namespace latl;

void test_block_ldlt()
{
    const int N = 3;
    const int n = 5;
    Matrix<N*n, N*n> bigm = Random();
    bigm = bigm.T() * bigm;
    Vector<N*n> bigv = Random();

    std::vector<Matrix<N,N> > bm(n*n);
    std::vector<Vector<N> > bv(n);

    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            bm[i*n+j] = slice<N,N>(bigm, i*N, j*N);
        }
        bv[i] = slice(i*N, N, bigv);
    }

    BlockMatrix<N,double> block_mat(bm.data(), n, n);

    BlockLDLT<N,double> ldlt;
    if (!ldlt.compute(block_mat)) {
        std::cerr << "error computing ldlt" << std::endl;
    }

    std::vector<Matrix<N,N> > bm_copy = bm;
    ldlt.reconstitute(block_mat);

    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            std::cerr << bm[i*n+j] - bm_copy[i*n+j] << std::endl;
        }
    }

    std::vector<Vector<N> > y(n);
    ldlt.inverse_times(bv, y);

    for (int i=0; i<n; ++i) {
        Vector<N> Ay(0.0);
        for (int j=0; j<n; ++j) {
            Ay += bm_copy[i*n+j] * y[j];
        }
        std::cerr << bv[i] - Ay << std::endl;
    }
}
#endif
