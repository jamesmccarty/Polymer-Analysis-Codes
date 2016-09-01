//  ALGORITHM COMPUTES DOMAIN RESOLVED CENTER-OF-MASS/CENTER-OF-MASS
//  TOTAL CORRELATION FUNCTION FOR A BINARY POLYMER MIXTURE
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
    char comd_[MXSZ], pcrd_[MXSZ], pfrm_[MXSZ],
        rf11_[MXSZ], rf12_[MXSZ], rf22_[MXSZ], rfTT_[MXSZ];
    double **amx, **amy, **amz, **bmx, **bmy, **bmz,
        *acx, *acy, *acz, *bcx, *bcy, *bcz, *h11, *h12, *h22, *hTT, ibox, idel,
        rmax, dist, xdis, ydis, zdis, delr, boxs, boxh, volu, rato, frx1;
    long i, m, n, k, t, bnmx, ste1, ste2, ply1, ply2, poly, nfrm;
    time_t t0, t1;
    FILE *fram, *cord, *f_11, *f_12, *f_22, *f_TT;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_HHPR50_ETHY50_T453_96/";     //  INPUT: DATA DIRECTORY NAME
    char fram_[] = "fram_HHPR50_ETHY50_T453_96";      //  INPUT: TIME FRAMES
    char cord_[] = "cord_HHPR50_ETHY50_T453_96";      //  INPUT: SITE COORDINATES
    char dirn_[] = "FBLN_HHPR50_ETHY50_T453_96/";     //  OUTPUT: RESULTS DIRECTORY NAME
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
    sprintf(rf11_, "%s%s", dirn_, f_11_);
    sprintf(rf12_, "%s%s", dirn_, f_12_);
    sprintf(rf22_, "%s%s", dirn_, f_22_);
    sprintf(rfTT_, "%s%s", dirn_, f_TT_);
    sprintf(comd_, "%s%s", "mkdir ", dirn_);
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
            //  READ MONOMER COORDINATES OF TYPE 1 CHAINS
            for(i=1; i<=ste1; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                amx[m][i] = xdis;
                amy[m][i] = ydis;
                amz[m][i] = zdis;
            }

            //  CALCULATE CENTER OF MASS COORDINATES OF TYPE 1 CHAINS
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(i=1; i<=ste1; i++)
            {
                xdis = xdis + amx[m][i];
                ydis = ydis + amy[m][i];
                zdis = zdis + amz[m][i];
            }
            acx[m] = xdis/(double)ste1;
            acy[m] = ydis/(double)ste1;
            acz[m] = zdis/(double)ste1;
        }

        //  PROCESS DATA OF TYPE 2 CHAINS
        for(m=1; m<=ply2; m++)
        {
            //  READ MONOMER COORDINATES OF TYPE 2 CHAINS
            for(i=1; i<=ste2; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                bmx[m][i] = xdis;
                bmy[m][i] = ydis;
                bmz[m][i] = zdis;
            }

            //  CALCULATE CENTER OF MASS COORDINATES OF TYPE 2 CHAINS
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(i=1; i<=ste2; i++)
            {
                xdis = xdis + bmx[m][i];
                ydis = ydis + bmy[m][i];
                zdis = zdis + bmz[m][i];
            }
            bcx[m] = xdis/(double)ste2;
            bcy[m] = ydis/(double)ste2;
            bcz[m] = zdis/(double)ste2;
        }

        //  CALCULATE TOTAL CORRELATION FUNCTIONS
        for(m=1; m<=ply1; m++)
        {
            for(n=1; n<=ply1; n++)
            {
                if(m!=n)
                {
                    //  CALCULATE hrc1c1
                    xdis = acx[m] - acx[n];
                    ydis = acy[m] - acy[n];
                    zdis = acz[m] - acz[n];
                    xdis = xdis - boxs*(double)nrst_long(xdis*ibox);
                    ydis = ydis - boxs*(double)nrst_long(ydis*ibox);
                    zdis = zdis - boxs*(double)nrst_long(zdis*ibox);
                    dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                    if(dist<rmax)
                    {
                        i = (long)(dist*idel);
                        h11[i] = h11[i] + 1.0;
                        hTT[i] = hTT[i] + 1.0;
                    }
                }
            }
            for(n=1; n<=ply2; n++)
            {
                //  CALCULATE hrc1c2
                xdis = acx[m] - bcx[n];
                ydis = acy[m] - bcy[n];
                zdis = acz[m] - bcz[n];
                xdis = xdis - boxs*(double)nrst_long(xdis*ibox);
                ydis = ydis - boxs*(double)nrst_long(ydis*ibox);
                zdis = zdis - boxs*(double)nrst_long(zdis*ibox);
                dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                if(dist<rmax)
                {
                    i = (long)(dist*idel);
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
                    xdis = bcx[m] - bcx[n];
                    ydis = bcy[m] - bcy[n];
                    zdis = bcz[m] - bcz[n];
                    xdis = xdis - boxs*(double)nrst_long(xdis*ibox);
                    ydis = ydis - boxs*(double)nrst_long(ydis*ibox);
                    zdis = zdis - boxs*(double)nrst_long(zdis*ibox);
                    dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                    if(dist<rmax)
                    {
                        i = (long)(dist*idel);
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
