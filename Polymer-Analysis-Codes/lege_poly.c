//  ROUTINE COMPUTES LEGENDRE POLYNOMIALS FOR BOND AND END-TO-END
//  VECTOR CORRELATORS.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutility.h"
#define MXSZ 500
void file_chck(FILE **, char *);
struct srcv
{
    double x, y, z;
};
int main()
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ],
        pp1b_[MXSZ], pp2b_[MXSZ], pp1e_[MXSZ], pp2e_[MXSZ],
        pa1b_[MXSZ], pa2b_[MXSZ], ps1b_[MXSZ], ps2b_[MXSZ];
    double **xcr1, **ycr1, **zcr1, **xcr2, **ycr2, **zcr2,
        **p1bm, **p2bm, *p1ev, *p2ev, *a1bv, *a2bv,
        lx_1, ly_1, lz_1, lx_2, ly_2, lz_2, dis1, dis2, dotp, tfac;
    long i, m, t1, t2, mono, poly, nfrm, lsze, fsze, *stmp;
    FILE *fram, *cord, *p1bd, *p2bd, *p1et, *p2et, *a1bd, *a2bd;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T450_0100/";            //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T450_0100";             //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "bmon_ETHY_T450_0100";             //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "LEGP_ETHY_T450_0100/";            //  OUTPUT: RESULTS DIRECTORY
    char dir1_[] = "LEG1_ETHY_T450_0100/";            //  OUTPUT: DIRECTORY RBND P1-FUNCTION
    char dir2_[] = "LEG2_ETHY_T450_0100/";            //  OUTPUT: DIRECTORY RBND P2-FUNCTION
    char p1bd_[] = "p1bd_ETHY_T450_0100";             //  OUTPUT: P1-FUNCTION, RBND
    char p2bd_[] = "p2bd_ETHY_T450_0100";             //  OUTPUT: P2-FUNCTION, RBND
    char p1et_[] = "p1et_ETHY_T450_0100";             //  OUTPUT: P1-FUNCTION, RETE
    char p2et_[] = "p2et_ETHY_T450_0100";             //  OUTPUT: P2-FUNCTION, RETE
    char a1bd_[] = "a1bd_ETHY_T450_0100";             //  OUTPUT: P1-FUNCTION, MEAN RBND
    char a2bd_[] = "a2bd_ETHY_T450_0100";             //  OUTPUT: P2-FUNCTION, MEAN RBND
    poly = 48;                                        //  NUMBER OF CHAINS
    mono = 100;                                       //  NUMBER OF MONOMERS PER CHAIN
    tfac = 0.01;                                      //  TIME CONVERSION FACTOR

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pp1e_, "%s%s", dirn_, p1et_);
    sprintf(pp2e_, "%s%s", dirn_, p2et_);
    sprintf(pa1b_, "%s%s", dirn_, a1bd_);
    sprintf(pa2b_, "%s%s", dirn_, a2bd_);
    sprintf(pp1b_, "%s%s%s", dirn_, dir1_, p1bd_);
    sprintf(pp2b_, "%s%s%s", dirn_, dir2_, p2bd_);
    sprintf(comd_, "mkdir %s", dirn_);
    system(comd_);
    sprintf(comd_, "mkdir %s%s", dirn_, dir1_);
    system(comd_);
    sprintf(comd_, "mkdir %s%s", dirn_, dir2_);
    system(comd_);

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&cord, pcrd_);
    fprintf(stdout, "Status of input files is satisfactory.\n");

    //  COUNT SIMULATION FRAMES
    fprintf(stdout, "Frames are being counted.\n");
    fram = fopen(pfrm_, "r");
    nfrm = 0;
    while((fscanf(fram, "%li", &i))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);
    nfrm = nfrm;

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");
    stmp = lvectr(1, nfrm);
    xcr1 = dmatrx(1, poly, 1, mono);
    ycr1 = dmatrx(1, poly, 1, mono);
    zcr1 = dmatrx(1, poly, 1, mono);
    xcr2 = dmatrx(1, poly, 1, mono);
    ycr2 = dmatrx(1, poly, 1, mono);
    zcr2 = dmatrx(1, poly, 1, mono);
    p1bm = dmatrx(1, nfrm, 1, mono-1);
    p2bm = dmatrx(1, nfrm, 1, mono-1);
    p1ev = dvectr(1, nfrm);
    p2ev = dvectr(1, nfrm);
    a1bv = dvectr(1, nfrm);
    a2bv = dvectr(1, nfrm);
    for(m=1; m<=nfrm; m++)
    {
        p1ev[m] = 0.0;
        p2ev[m] = 0.0;
        a1bv[m] = 0.0;
        a2bv[m] = 0.0;
        for(i=1; i<mono; i++)
        {
            p1bm[m][i] = 0.0;
            p2bm[m][i] = 0.0;
        }
    }

    //  STORE SIMULATION FRAMES
    fprintf(stdout, "Frames are being stored.\n");
    fram = fopen(pfrm_, "r");
    for(i=1; i<=nfrm; i++)
    {
        fscanf(fram, "%li", &stmp[i]);
    }
    fclose(fram);

    //  COMPUTE LEGENDRE POLYNOMIALS
    fprintf(stdout, "Legendre polynomials are being calculated.\n");
    lsze = sizeof(struct srcv);
    fsze = mono*poly*lsze;
    cord = fopen(pcrd_, "r");
    for(t1=1; t1<=nfrm; t1++)
    {
        si = time(NULL);
        fprintf(stdout, "CYCLE %5li OF %5li\t", t1, nfrm);
        fseek(cord, (t1-1)*fsze, SEEK_SET);
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
                fread(&cont, lsze, 1, cord);
                xcr1[m][i] = cont.x;
                ycr1[m][i] = cont.y;
                zcr1[m][i] = cont.z;
            }
        }
        for(t2=t1; t2<=nfrm; t2++)
        {
            fseek(cord, (t2-1)*fsze, SEEK_SET);
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<=mono; i++)
                {
                    fread(&cont, lsze, 1, cord);
                    xcr2[m][i] = cont.x;
                    ycr2[m][i] = cont.y;
                    zcr2[m][i] = cont.z;
                }
            }
            for(m=1; m<=poly; m++)
            {
                lx_1 = xcr1[m][mono] - xcr1[m][1];
                ly_1 = ycr1[m][mono] - ycr1[m][1];
                lz_1 = zcr1[m][mono] - zcr1[m][1];
                lx_2 = xcr2[m][mono] - xcr2[m][1];
                ly_2 = ycr2[m][mono] - ycr2[m][1];
                lz_2 = zcr2[m][mono] - zcr2[m][1];
                dis1 = sqrt(lx_1*lx_1 + ly_1*ly_1 + lz_1*lz_1);
                dis2 = sqrt(lx_2*lx_2 + ly_2*ly_2 + lz_2*lz_2);
                dotp = (lx_1*lx_2 + ly_1*ly_2 + lz_1*lz_2)/(dis1*dis2);
                p1ev[t2-t1+1] = p1ev[t2-t1+1] + dotp;
                p2ev[t2-t1+1] = p2ev[t2-t1+1] + dotp*dotp;
                for(i=1; i<mono; i++)
                {
                    lx_1 = xcr1[m][i+1] - xcr1[m][i];
                    ly_1 = ycr1[m][i+1] - ycr1[m][i];
                    lz_1 = zcr1[m][i+1] - zcr1[m][i];
                    lx_2 = xcr2[m][i+1] - xcr2[m][i];
                    ly_2 = ycr2[m][i+1] - ycr2[m][i];
                    lz_2 = zcr2[m][i+1] - zcr2[m][i];
                    dis1 = sqrt(lx_1*lx_1 + ly_1*ly_1 + lz_1*lz_1);
                    dis2 = sqrt(lx_2*lx_2 + ly_2*ly_2 + lz_2*lz_2);
                    dotp = (lx_1*lx_2 + ly_1*ly_2 + lz_1*lz_2)/(dis1*dis2);
                    p1bm[t2-t1+1][i] = p1bm[t2-t1+1][i] + dotp;
                    p2bm[t2-t1+1][i] = p2bm[t2-t1+1][i] + dotp*dotp;
                    a1bv[t2-t1+1] = a1bv[t2-t1+1] + dotp;
                    a2bv[t2-t1+1] = a2bv[t2-t1+1] + dotp*dotp;
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }

    //  NORMALIZE AND PRINT RESULTS FOR LEGENDRE POLYNOMIALS
    fprintf(stdout, "Results are being printed.\n");
    a1bd = fopen(pa1b_, "w");
    a2bd = fopen(pa2b_, "w");
    p1et = fopen(pp1e_, "w");
    p2et = fopen(pp2e_, "w");
    for(i=1; i<=nfrm; i++)
    {
        p1ev[i] = p1ev[i]/(double)((nfrm-i+1)*poly);
        p2ev[i] = 1.5*p2ev[i]/(double)((nfrm-i+1)*poly)-0.5;
        a1bv[i] = a1bv[i]/(double)((nfrm-i+1)*(mono-1)*poly);
        a2bv[i] = 1.5*a2bv[i]/(double)((nfrm-i+1)*(mono-1)*poly)-0.5;
        fprintf(p1et, "% 20.16E\t% 20.16E\n", tfac*stmp[i], p1ev[i]);
        fprintf(p2et, "% 20.16E\t% 20.16E\n", tfac*stmp[i], p2ev[i]);
        fprintf(a1bd, "% 20.16E\t% 20.16E\n", tfac*stmp[i], a1bv[i]);
        fprintf(a2bd, "% 20.16E\t% 20.16E\n", tfac*stmp[i], a2bv[i]);
    }
    fclose(p1et);
    fclose(p2et);
    fclose(a1bd);
    fclose(a2bd);
    for(m=1; m<mono; m++)
    {
        sprintf(ps1b_, "%s_B%03li", pp1b_, m);
        sprintf(ps2b_, "%s_B%03li", pp2b_, m);
        p1bd = fopen(ps1b_, "w");
        p2bd = fopen(ps2b_, "w");
        for(i=1; i<=nfrm; i++)
        {
            p1bm[i][m] = p1bm[i][m]/(double)((nfrm-i+1)*poly);
            p2bm[i][m] = 1.5*p2bm[i][m]/(double)((nfrm-i+1)*poly)-0.5;
            fprintf(p1bd, "% 20.16E\t% 20.16E\n", tfac*stmp[i], p1bm[i][m]);
            fprintf(p2bd, "% 20.16E\t% 20.16E\n", tfac*stmp[i], p2bm[i][m]);
        }
        fclose(p1bd);
        fclose(p2bd);
    }

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
    return(0);
}

//  FILE EXISTENCE CHECKER
void file_chck(FILE **fptr, char fstr_[])
{
    if(!(*fptr = fopen(fstr_, "r")))
    {
        fprintf(stdout, "Failed to open file %s\n", fstr_);
        exit(1);
    }
    else
    {
        fclose(*fptr);
    }
}
