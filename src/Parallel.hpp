//
//  Parallel.hpp
//  Muffin
//
//  Threads of MUFFIN, set by the environment variable MUFFIN_THREADS (1 by default: the serial code, without
//  OpenMP region). With more threads the loops over the cells and the faces are shared among them; every cell is
//  computed by the same operations and the sums are made in the original order, so the results are the same to
//  the bit whatever the number of threads.
//

#ifndef Parallel_hpp
#define Parallel_hpp

#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

inline int muffinThreads()
{
#ifdef _OPENMP
    static const int n = [] {
        const char* s = std::getenv("MUFFIN_THREADS");
        const int v = s ? std::atoi(s) : 1;
        return v > 1 ? v : 1;
    }();
    return n;
#else
    return 1;
#endif
}

inline int muffinThreadId()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// f(begin, end, thread) on contiguous parts of [begin, end), one per thread; with one thread, f(begin, end, 0).
template <class F>
inline void parallelRange(int begin, int end, F&& f)
{
#ifdef _OPENMP
    const int nt = muffinThreads();
    if (nt > 1) {
        #pragma omp parallel num_threads(nt)
        {
            const int t = omp_get_thread_num(), n = omp_get_num_threads();
            const long len = end - begin;
            const int b = begin + static_cast<int>(len*t/n), e = begin + static_cast<int>(len*(t + 1)/n);
            if (b < e) f(b, e, t);
        }
        return;
    }
#endif
    if (begin < end) f(begin, end, 0);
}

#endif /* Parallel_hpp */
