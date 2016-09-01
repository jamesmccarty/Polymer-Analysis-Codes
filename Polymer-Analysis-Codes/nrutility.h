#ifndef _NRUTILITY_H_
#define _NRUTILITY_H_
#endif

long   *lvectr(long, long);
long  **lmatrx(long, long, long, long);
long ***ltensr(long, long, long, long, long, long);

void free_lvectr(long   *, long, long);
void free_lmatrx(long  **, long, long, long, long);
void free_ltensr(long ***, long, long, long, long, long, long);

double   *dvectr(long, long);
double  **dmatrx(long, long, long, long);
double ***dtensr(long, long, long, long, long, long);

void free_dvectr(double   *, long, long);
void free_dmatrx(double  **, long, long, long, long);
void free_dtensr(double ***, long, long, long, long, long, long);
