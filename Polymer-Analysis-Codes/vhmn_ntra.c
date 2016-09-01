//  ROUTINE COMPUTES INTRAMOLECULAR MONOMER VAN HOVE FUNCTION
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutility.h"
#define pi M_PI
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
    char comd_[MXSZ], pgmr_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], pstp_[MXSZ];
    double **xcr1, **ycr1, **zcr1, **xcr2, **ycr2, **zcr2, **hist,
        xdis, ydis, zdis, dist, boxs, delr, invd, invb,
        rmax, shex, shin, volu, scal;
    long *delt, i, j, k, m, t1, t2, ibin, mono, poly, nbin,
        nfrm, cnts, lsze, fsze, tnum, *stmp;
    FILE *cord, *fram, *gmrt, *step;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T450_0100/";            //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T450_0100";             //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "bmon_ETHY_T450_0100";             //  INPUT: BINARY SITE COORDINATES
    char fstp_[] = "fstp_ETHY_T450_0100";             //  INPUT: TIME INTERVALS
    char dirn_[] = "VHRT_ETHY_T450_0100/";            //  OUTPUT: RESULTS DIRECTORY
    char gmrt_[] = "gmra_ETHY_T450_0100";             //  OUTPUT: VAN HOVE FUNCTION
    poly = 48;                                        //  NUMBER OF CHAINS
    mono = 100;                                       //  NUMBER OF MONOMERS PER CHAIN
    delr = 0.04;                                      //  HISTOGRAM BIN WIDTH
    boxs = 52.9;                                      //  MEAN BOX DIMENSION
    scal = 2.0;                                       //  BOX DIMENSION SCALING FACTOR

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pstp_, "%s%s", base_, fstp_);
    sprintf(comd_, "mkdir %s", dirn_);
    system(comd_);

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&cord, pcrd_);
    file_chck(&step, pstp_);
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

    //  COUNT TIME INTERVALS
    fprintf(stdout, "Intervals are being counted.\n");
    step = fopen(pstp_, "r");
    cnts = 0;
    while((fscanf(step, "%li", &i))!=EOF)
    {
        cnts++;
    }
    fclose(step);
    cnts = cnts;

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");
    nbin = (long)(boxs*scal/delr);
    stmp = lvectr(1, nfrm);
    delt = lvectr(1, cnts);
    hist = dmatrx(1, cnts, 0, nbin-1);
    xcr1 = dmatrx(1, poly, 1, mono);
    ycr1 = dmatrx(1, poly, 1, mono);
    zcr1 = dmatrx(1, poly, 1, mono);
    xcr2 = dmatrx(1, poly, 1, mono);
    ycr2 = dmatrx(1, poly, 1, mono);
    zcr2 = dmatrx(1, poly, 1, mono);
    for(m=1; m<=cnts; m++)
    {
        for(i=0; i<nbin; i++)
        {
            hist[m][i] = 0.0;
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

    //  STORE TIME INTERVALS
    fprintf(stdout, "Time intervals are being stored.\n");
    step = fopen(pstp_, "r");
    for(i=1; i<=cnts; i++)
    {
        fscanf(step, "%li", &delt[i]);
    }
    fclose(step);

    //  COMPUTE HISTOGRAM FOR VAN HOVE FUNCTION
    fprintf(stdout, "Histogram is being calculated.\n");
    invd = 1.0/delr;
    invb = 1.0/boxs;
    rmax = (double)(nbin*delr);
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
            for(k=1; k<=cnts; k++)
            {
                if((stmp[t2]-stmp[t1])==delt[k])
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
                        for(i=1; i<=mono; i++)
                        {
                            for(j=1; j<=mono; j++)
                            {
                                if(i!=j)
                                {
                                    xdis = xcr2[m][j] - xcr1[m][i];
                                    ydis = ycr2[m][j] - ycr1[m][i];
                                    zdis = zcr2[m][j] - zcr1[m][i];
                                    dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                                    if(dist<rmax)
                                    {
                                        ibin = (long)(dist*invd);
                                        hist[k][ibin] = hist[k][ibin] + 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR VAN HOVE FUNCTION
    fprintf(stdout, "Results are being printed.\n");
    for(m=1; m<=cnts; m++)
    {
        sprintf(pgmr_, "%s%s_D%li", dirn_, gmrt_, delt[m]);
        gmrt = fopen(pgmr_, "w");
        tnum = 0;
        for(t1=1; t1<=nfrm; t1++)
        {
            for(t2=t1; t2<=nfrm; t2++)
            {
                if((stmp[t2]-stmp[t1])==delt[m])
                {
                    tnum++;
                }
            }
        }
        for(i=0; i<nbin; i++)
        {
            shex = (i+1)*(i+1)*(i+1);
            shin = i*i*i;
            volu = 4.0/3.0*pi*(shex-shin)*delr*delr*delr;
            hist[m][i] = hist[m][i]
                    / (double)(mono*poly*mono*poly*tnum*volu*invb*invb*invb);
            fprintf(gmrt, "% 20.16E\t% 20.16E\n", (double)((i+0.5)*delr), hist[m][i]);
        }
        fclose(gmrt);
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
