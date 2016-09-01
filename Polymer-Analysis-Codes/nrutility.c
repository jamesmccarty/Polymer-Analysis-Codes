//  LIBRARY FOR MEMORY ALLOCATION
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1

//  ALLOCATE VECTOR OF DOUBLES WITH RANGE v[nrl...nrh]
double *dvectr(long nrl, long nrh)
{
    long nrow;
    double *v;

    nrow = nrh - nrl + 1;
    v = (double *)malloc((size_t)((nrow+NR_END)*sizeof(double)));
    if(!v)
    {
        printf("Allocation failure in dvectr().\n");
        exit(1);
    }
    return(v-nrl+NR_END);
}

//  FREE VECTOR OF DOUBLES FROM MEMORY
void free_dvectr(double *v, long nrl, long nrh)
{
    free((char *)(v+nrl-NR_END));
}

//  ALLOCATE MATRIX OF DOUBLES WITH RANGE m[nrl...nrh][ncl...nch]
double **dmatrx(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow, ncol;
    double **m;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;

    //  ALLOCATE POINTERS TO ROWS
    m = (double **)malloc((size_t)((nrow+NR_END)*sizeof(double *)));
    if(!m)
    {
        printf("Allocation failure 1 in dmatrx().\n");
        exit(1);
    }
    m += NR_END;
    m -= nrl;

    //  ALLOCATE ROWS AND SET POINTERS TO THESE
    m[nrl] = (double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if(!m[nrl])
    {
        printf("Allocation failure 2 in dmatrx().\n");
        exit(1);
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++)
    {
        m[i] = m[i-1] + ncol;
    }

    //  RETURN POINTER TO ARRAY OF POINTERS TO ROWS
    return(m);
}

//  FREE MATRIX OF DOUBLES FROM MEMORY
void free_dmatrx(double **m, long nrl, long nrh, long ncl, long nch)
{
    free((char *)(m[nrl]+ncl-NR_END));
    free((char *)(m+nrl-NR_END));
}

//  ALLOCATE TENSOR OF DOUBLES WITH RANGE t[nrl...nrh][ncl...nch][ndl...ndh]
double ***dtensr(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    long i, j, nrow, ncol, ndim;
    double ***t;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;
    ndim = ndh - ndl + 1;

    //  ALLOCATE POINTERS TO POINTERS TO ROWS
    t = (double ***)malloc((size_t)((nrow+NR_END)*sizeof(double **)));
    if(!t)
    {
        printf("Allocation failure 1 in dtensr().\n");
        exit(1);
    }
    t += NR_END;
    t -= nrl;

    //  ALLOCATE POINTERS TO ROWS AND SET POINTERS TO THESE
    t[nrl] = (double **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double *)));
    if(!t[nrl])
    {
        printf("Allocation failure 2 in dtensr().\n");
        exit(1);
    }
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    //  ALLOCATE ROWS AND SET POINTERS TO THESE
    t[nrl][ncl] = (double *)malloc((size_t)((nrow*ncol*ndim+NR_END)*sizeof(double)));
    if(!t[nrl][ncl])
    {
        printf("Allocation failure 3 in dtensr().\n");
        exit(1);
    }

    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1; j<=nch; j++)
    {
        t[nrl][j] = t[nrl][j-1] + ndim;
    }
    for(i=nrl+1; i<=nrh; i++)
    {
        t[i] = t[i-1] + ncol;
        t[i][ncl] = t[i-1][ncl] + ncol*ndim;
        for(j=ncl+1; j<=nch; j++)
        {
            t[i][j] = t[i][j-1] + ndim;
        }
    }

    //  RETURN POINTER TO ARRAY OF POINTERS TO ROWS
    return(t);
}

//  FREE TENSOR OF DOUBLES FROM MEMORY
void free_dtensr(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    free((char *)(t[nrl][ncl]+ndl-NR_END));
    free((char *)(t[nrl]+ncl-NR_END));
    free((char *)(t+nrl-NR_END));
}

//  ALLOCATE VECTOR OF LONGS WITH RANGE v[nrl...nrh]
long *lvectr(long nrl, long nrh)
{
    long nrow;
    long *v;

    nrow = nrh - nrl + 1;
    v = (long *)malloc((size_t)((nrow+NR_END)*sizeof(long)));
    if(!v)
    {
        printf("Allocation failure in lvectr().\n");
        exit(1);
    }
    return(v-nrl+NR_END);
}

//  FREE VECTOR OF LONGS FROM MEMORY
void free_lvectr(long *v, long nrl, long nrh)
{
    free((char *)(v+nrl-NR_END));
}

//  ALLOCATE MATRIX OF LONGS WITH RANGE m[nrl...nrh][ncl...nch]
long **lmatrx(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow, ncol;
    long **m;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;

    //  ALLOCATE POINTERS TO ROWS
    m = (long **)malloc((size_t)((nrow+NR_END)*sizeof(long *)));
    if(!m)
    {
        printf("Allocation failure 1 in lmatrx().\n");
        exit(1);
    }
    m += NR_END;
    m -= nrl;

    //  ALLOCATE ROWS AND SET POINTERS TO THESE
    m[nrl] = (long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
    if(!m[nrl])
    {
        printf("Allocation failure 2 in lmatrx().\n");
        exit(1);
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++)
    {
        m[i] = m[i-1] + ncol;
    }

    //  RETURN POINTER TO ARRAY OF POINTERS TO ROWS
    return(m);
}

//  FREE MATRIX OF LONGS FROM MEMORY
void free_lmatrx(long **m, long nrl, long nrh, long ncl, long nch)
{
    free((char *)(m[nrl]+ncl-NR_END));
    free((char *)(m+nrl-NR_END));
}

//  ALLOCATE TENSOR LONGS WITH RANGE t[nrl...nrh][ncl...nch][ndl...ndh]
long ***ltensr(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    long i, j, nrow, ncol, ndim;
    long ***t;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;
    ndim = ndh - ndl + 1;

    //  ALLOCATE POINTERS TO POINTERS TO ROWS
    t = (long ***)malloc((size_t)((nrow+NR_END)*sizeof(long **)));
    if(!t)
    {
        printf("Allocation failure 1 in ltensr().\n");
        exit(1);
    }
    t += NR_END;
    t -= nrl;

    //  ALLOCATE POINTERS TO ROWS AND SET POINTERS TO THESE
    t[nrl] = (long **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long *)));
    if(!t[nrl])
    {
        printf("Allocation failure 2 in ltensr().\n");
        exit(1);
    }
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    //  ALLOCATE ROWS AND SET POINTERS TO THESE
    t[nrl][ncl] = (long *)malloc((size_t)((nrow*ncol*ndim+NR_END)*sizeof(long)));
    if(!t[nrl][ncl])
    {
        printf("Allocation failure 3 in ltensr().\n");
        exit(1);
    }

    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1; j<=nch; j++)
    {
        t[nrl][j] = t[nrl][j-1] + ndim;
    }
    for(i=nrl+1; i<=nrh; i++)
    {
        t[i] = t[i-1] + ncol;
        t[i][ncl] = t[i-1][ncl] + ncol*ndim;
        for(j=ncl+1; j<=nch; j++)
        {
            t[i][j] = t[i][j-1] + ndim;
        }
    }

    //  RETURN POINTER TO ARRAY OF POINTERS TO ROWS
    return(t);
}

//  FREE TENSOR OF LONGS FROM MEMORY
void free_ltensr(long ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    free((char *)(t[nrl][ncl]+ndl-NR_END));
    free((char *)(t[nrl]+ncl-NR_END));
    free((char *)(t+nrl-NR_END));
}
