//  ALGORITHM COMPUTES DOMAIN RESOLVED CENTER OF MASS-CENTER OF MASS
//  TOTAL CORRELATION FUNCTION (REAL SPACE) FOR A BINARY POLYMER MIXUTRE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutility.h"
#define pi M_PI
#define MXSZ 500
void file_chck(FILE **, char *);
long nrst_long(double);
int main()
{
    double **amx, **amy, **amz, **bmx, **bmy, **bmz,
        *acx, *acy, *acz, *bcx, *bcy, *bcz, *h11, *h12, *h22, *hTT,
        corx, cory, corz, cnrx, cnry, cnrz, ibox, idel, rmax,
        crix, criy, criz, delr, boxs, boxh, volu, rato, frx1;
    char comd_[MXSZ], dirn_[MXSZ], pcrd_[MXSZ], pfrm_[MXSZ],
        rf11_[MXSZ], rf12_[MXSZ], rf22_[MXSZ], rfTT_[MXSZ];
    long i, m, n, k, t, bnmx, ste1, ste2, ply1, ply2, poly, nfrm;
    time_t t0, t1;
    FILE *fram, *cord, *f_11, *f_12, *f_22, *f_TT;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_HHPR50_ETHY50_T453_96/";     //  INPUT: DATA DIRECTORY NAME
    char fram_[] = "fram_HHPR50_ETHY50_T453_96";      //  INPUT: TIME FRAMES
    char cord_[] = "cord_HHPR50_ETHY50_T453_96";      //  INPUT: SITE COORDINATES
    char dirp_[] = "FBLN_HHPR50_ETHY50_T453_96/";     //  OUTPUT: RESULTS DIRECTORY NAME
    char f_11_[] = "hc1c1r_HHPR50_ETHY50_T453_96";    //  OUTPUT: Hc1c1R
    char f_12_[] = "hc1c2r_HHPR50_ETHY50_T453_96";    //  OUTPUT: Hc1c2R
    char f_22_[] = "hc2c2r_HHPR50_ETHY50_T453_96";    //  OUTPUT: Hc2c2R
    char f_TT_[] = "hcTcTr_HHPR50_ETHY50_T453_96";    //  OUTPUT: HcTcTR
    ste1 = 96;                                        //  TYPE 1 NUMBER OF SITES PER CHAIN
    ste2 = 96;                                        //  TYPE 2 NUMBER OF SITES PER CHAIN
    frx1 = 0.50;                                      //  TYPE 1 FRACTION OF CHAINS
    poly = 1600;                                      //  NUMBER OF POLYMERS
    delr = 0.01*sqrt((ste1*ste1+ste2*ste2)/2);        //  BIN WIDTH
    boxs = 166.612;                                   //  LENGTH OF SIMULATION CELL

    //  AUXILIARY PARAMETERS
    ibox = 1.0/boxs;                                  //  INVERSE BOX LENGTH
    idel = 1.0/delr;                                  //  INVERSE BIN WIDTH
    boxh = boxs/2.0;                                  //  HALVED LINER DIMENSION OF SIMULATION BOX
    bnmx = (long)(boxh/delr);                         //  MAXIMUM NUMBER OF BINS
    rmax = (double)(bnmx*delr);                       //  MAXIMUM RADIAL DISTANCE
    volu = 4.0/3.0*pi*pow(delr/boxs, 3.0);            //  ENSEMBLE VOLUME FACTOR
    ply1 = (long)(poly*frx1);                         //  TYPE 1 NUMBER OF CHAINS
    ply2 = (long)(poly*(1.0-frx1));                   //  TYPE 2 NUMBER OF CHAINS

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(rf11_, "%s%s", dirp_, f_11_);
    sprintf(rf12_, "%s%s", dirp_, f_12_);
    sprintf(rf22_, "%s%s", dirp_, f_22_);
    sprintf(rfTT_, "%s%s", dirp_, f_TT_);
    sprintf(comd_, "%s%s", "mkdir ", dirp_);
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
    fprintf(stdout, "Memory resources are begin allocated.\n");
    amx = dmatrx(1,ply1,1,ste1);
    amy = dmatrx(1,ply1,1,ste1);
    amz = dmatrx(1,ply1,1,ste1);
    bmx = dmatrx(1,ply2,1,ste2);
    bmy = dmatrx(1,ply2,1,ste2);
    bmz = dmatrx(1,ply2,1,ste2);
    acx = dvectr(1,ply1);
    acy = dvectr(1,ply1);
    acz = dvectr(1,ply1);
    bcx = dvectr(1,ply2);
    bcy = dvectr(1,ply2);
    bcz = dvectr(1,ply2);
    h11 = dvectr(0,bnmx);
    h12 = dvectr(0,bnmx);
    h22 = dvectr(0,bnmx);
    hTT = dvectr(0,bnmx);
    for(i=0; i<=bnmx; i++)
    {
        h11[i] = 0.0;
        h12[i] = 0.0;
        h22[i] = 0.0;
        hTT[i] = 0.0;
    }

    //  CALCULATE CORRELATION FUNCTIONS
    fprintf(stdout, "Correlation functions are being calculated.\n");
    cord = fopen(pcrd_, "r");
    for(t=1; t<=nfrm; t++)
    {
        fprintf(stdout, "FRAME %5li OF %5li\t", t, nfrm);
        t0 = time(NULL);

        //  PROCESS DATA OF TYPE 1 CHAINS
        for(m=1; m<=ply1; m++)
        {
            corx = 0.0;
            cory = 0.0;
            corz = 0.0;

            //  READ MONOMER COORDINATES OF TYPE 1 CHAINS
            for(i=1; i<=ste1; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &crix, &criy, &criz);
                amx[m][i] = crix;
                amy[m][i] = criy;
                amz[m][i] = criz;
                corx = corx + crix;
                cory = cory + criy;
                corz = corz + criz;
            }
            corx = corx/(double)ste1;
            cory = cory/(double)ste1;
            corz = corz/(double)ste1;

            //  APPLY SHIFT CORRECTIONS
            cnrx = corx - boxs*nrst_long(corx*ibox);
            cnry = cory - boxs*nrst_long(cory*ibox);
            cnrz = corz - boxs*nrst_long(corz*ibox);
            for(i=1; i<=ste1; i++)
            {
                amx[m][i] = amx[m][i] - corx + cnrx;
                amy[m][i] = amy[m][i] - cory + cnry;
                amz[m][i] = amz[m][i] - corz + cnrz;
            }

            //  CALCULATE CENTER OF MASS COORDINATES OF TYPE 1 CHAINS
            crix = 0.0;
            criy = 0.0;
            criz = 0.0;
            for(i=1; i<=ste1; i++)
            {
                crix = crix + amx[m][i];
                criy = criy + amy[m][i];
                criz = criz + amz[m][i];
            }
            acx[m] = crix/(double)ste1;
            acy[m] = criy/(double)ste1;
            acz[m] = criz/(double)ste1;
        }

        //  PROCESS DATA OF TYPE 2 CHAINS
        for(m=1; m<=ply2; m++)
        {
            corx = 0.0;
            cory = 0.0;
            corz = 0.0;

            //  READ MONOMER COORDINATES OF TYPE 2 CHAINS
            for(i=1; i<=ste2; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &crix, &criy, &criz);
                bmx[m][i] = crix;
                bmy[m][i] = criy;
                bmz[m][i] = criz;
                corx = corx + crix;
                cory = cory + criy;
                corz = corz + criz;
            }
            corx = corx/(double)ste2;
            cory = cory/(double)ste2;
            corz = corz/(double)ste2;

            //  APPLY SHIFT CORRECTIONS
            cnrx = corx - boxs*nrst_long(corx*ibox);
            cnry = cory - boxs*nrst_long(cory*ibox);
            cnrz = corz - boxs*nrst_long(corz*ibox);
            for(i=1; i<=ste2; i++)
            {
                bmx[m][i] = bmx[m][i] - corx + cnrx;
                bmy[m][i] = bmy[m][i] - cory + cnry;
                bmz[m][i] = bmz[m][i] - corz + cnrz;
            }

            //  CALCULATE CENTER OF MASS COORDINATES OF TYPE 2 CHAINS
            crix = 0.0;
            criy = 0.0;
            criz = 0.0;
            for(i=1; i<=ste2; i++)
            {
                crix = crix + bmx[m][i];
                criy = criy + bmy[m][i];
                criz = criz + bmz[m][i];
            }
            bcx[m] = crix/(double)ste2;
            bcy[m] = criy/(double)ste2;
            bcz[m] = criz/(double)ste2;
        }

        //  CALCULATE TOTAL CORRELATION FUNCTIONS
        for(m=1; m<=ply1; m++)
        {
            for(n=1; n<=ply1; n++)
            {
                if(m!=n)
                {
                    //  CALCULATE hrc1c1
                    crix = acx[m] - acx[n];
                    criy = acy[m] - acy[n];
                    criz = acz[m] - acz[n];
                    crix = crix - boxs*nrst_long(crix*ibox);
                    criy = criy - boxs*nrst_long(criy*ibox);
                    criz = criz - boxs*nrst_long(criz*ibox);
                    crix = sqrt(crix*crix + criy*criy + criz*criz);
                    if(crix<rmax)
                    {
                        i = (long)(crix*idel);
                        h11[i] = h11[i] + 1.0;
                        hTT[i] = hTT[i] + 1.0;
                    }
                }
            }
            for(n=1; n<=ply2; n++)
            {
                //  CALCULATE hrc1c2
                crix = acx[m] - bcx[n];
                criy = acy[m] - bcy[n];
                criz = acz[m] - bcz[n];
                crix = crix - boxs*nrst_long(crix*ibox);
                criy = criy - boxs*nrst_long(criy*ibox);
                criz = criz - boxs*nrst_long(criz*ibox);
                crix = sqrt(crix*crix + criy*criy + criz*criz);
                if(crix<rmax)
                {
                    i = (long)(crix*idel);
                    h12[i] = h12[i] + 1.0;
                    hTT[i] = hTT[i] + 2.0;
                }
            }
        }
        for(m=1; m<=ply2; m++)
        {
            for(n=1; n<=ply2; n++)
            {
                if(m!=n)
                {
                    //  CALCULATE hrc2c2
                    crix = bcx[m] - bcx[n];
                    criy = bcy[m] - bcy[n];
                    criz = bcz[m] - bcz[n];
                    crix = crix - boxs*nrst_long(crix*ibox);
                    criy = criy - boxs*nrst_long(criy*ibox);
                    criz = criz - boxs*nrst_long(criz*ibox);
                    crix = sqrt(crix*crix + criy*criy + criz*criz);
                    if(crix<rmax)
                    {
                        i = (long)(crix*idel);
                        h22[i] = h22[i] + 1.0;
                        hTT[i] = hTT[i] + 1.0;
                    }
                }
            }
        }
        t1 = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", t1-t0);
    }
    fclose(cord);

    //  PRINT FORM FACTOR RESULTS
    fprintf(stdout, "Results are being printed.\n");
    f_11 = fopen(rf11_, "w");
    f_12 = fopen(rf12_, "w");
    f_22 = fopen(rf22_, "w");
    f_TT = fopen(rfTT_, "w");
    for(i=0; i<bnmx; i++)
    {
        rato = volu*(double)((i+1)*(i+1)*(i+1)-i*i*i);
        h11[i] = h11[i]/rato/(double)(ply1*ply1*nfrm);
        h12[i] = h12[i]/rato/(double)(ply1*ply2*nfrm);
        h22[i] = h22[i]/rato/(double)(ply2*ply2*nfrm);
        hTT[i] = hTT[i]/rato/(double)(poly*poly*nfrm);
        fprintf(f_11, "% 20.16E\t% 20.16E\n", (i+0.5)*delr, h11[i]-1.0);
        fprintf(f_12, "% 20.16E\t% 20.16E\n", (i+0.5)*delr, h12[i]-1.0);
        fprintf(f_22, "% 20.16E\t% 20.16E\n", (i+0.5)*delr, h22[i]-1.0);
        fprintf(f_TT, "% 20.16E\t% 20.16E\n", (i+0.5)*delr, hTT[i]-1.0);
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

//  COMPUTE NEAREST INTEGER
long nrst_long(double argm)
{
    if(argm>0)
    {
        return (long)(argm + 0.5);
    }
    else
    {
        return (long)(argm - 0.5);
    }
}
