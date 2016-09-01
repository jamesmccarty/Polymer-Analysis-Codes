//  ROUTINE COMPUTES MOLECULAR LENGTH STATISTICS
//  FOR A BINARY POLYMER MIXTURE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutility.h"
#define MXSZ 500
void file_chck(FILE **, char *);
int main()
{
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ],
        prg1_[MXSZ], prb1_[MXSZ], pre1_[MXSZ],
        prg2_[MXSZ], prb2_[MXSZ], pre2_[MXSZ];
    double **amx, **amy, **amz, **bmx, **bmy, **bmz,
        **xam, **yam, **zam, **xbm, **ybm, **zbm,
        *acx, *acy, *acz, *bcx, *bcy, *bcz,
        xdis, ydis, zdis, grba, grbb, grga, grgb, grea, greb, dist, frx1;
    long a, i, j, k, nfrm, step, poly,
        ply1, mon1, bds1, ste1, ply2, mon2, bds2, ste2;
    FILE *cord, *fram, *frb1, *frb2, *frg1, *frg2, *fre1, *fre2;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_HHPR50_ISBU50_T453_96/";     //  INPUT: DATA DIRECTORY
    char cord_[] = "cord_HHPR50_ISBU50_T453_96";      //  INPUT: SITE COORDINATES
    char fram_[] = "fram_HHPR50_ISBU50_T453_96";      //  INPUT: SIMULATION TIME STEPS
    char dirn_[] = "LENG_HHPR50_ISBU50_T453_96/";     //  OUTPUT: RESULTS DIRECTORY
    char fng1_[] = "rgyr_HHPR50_T453_96";             //  OUTPUT: TYPE 1 MOLECULAR GYRATION RADIUS
    char fnb1_[] = "rbnd_HHPR50_T453_96";             //  OUTPUT: TYPE 1 MONOMER BOND LENGTH
    char fne1_[] = "rext_HHPR50_T453_96";             //  OUTPUT: TYPE 1 POLYMER EXTENSION
    char fng2_[] = "rgyr_ISBU50_T453_96";             //  OUTPUT: TYPE 2 MOLECULAR GYRATION RADIUS
    char fnb2_[] = "rbnd_ISBU50_T453_96";             //  OUTPUT: TYPE 2 MONOMER BOND LENGTH
    char fne2_[] = "rext_ISBU50_T453_96";             //  OUTPUT: TYPE 2 POLYMER EXTENSION
    ste1 = 96;                                        //  TYPE 1 NUMBER OF SITES PER CHAIN
    ste2 = 96;                                        //  TYPE 2 NUMBER OF SITES PER CHAIN
    bds1 = 6;                                         //  TYPE 1 NUMBER OF SITES IN MONOMER
    bds2 = 4;                                         //  TYPE 2 NUMBER OF SITES IN MONOMER
    frx1 = 0.50;                                      //  TYPE 1 FRACTION OF CHAINS
    poly = 1600;                                      //  NUMBER OF POLYMERS

    //  AUXILIARY PARAMETERS
    ply1 = (long)(poly*frx1);                         //  TYPE 1 NUMBER OF CHAINS
    mon1 = ste1/bds1;                                 //  TYPE 1 NUMBER OF MONOMERS IN CHAIN
    ply2 = (long)(poly*(1.0-frx1));                   //  TYPE 2 NUMBER OF CHAINS
    mon2 = ste2/bds2;                                 //  TYPE 2 NUMBER OF MONOMERS IN CHAIN

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(prg1_, "%s%s", dirn_, fng1_);
    sprintf(prb1_, "%s%s", dirn_, fnb1_);
    sprintf(pre1_, "%s%s", dirn_, fne1_);
    sprintf(prg2_, "%s%s", dirn_, fng2_);
    sprintf(prb2_, "%s%s", dirn_, fnb2_);
    sprintf(pre2_, "%s%s", dirn_, fne2_);
    sprintf(comd_, "mkdir %s", dirn_);
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
    while((fscanf(fram, "%li", &step))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);
    nfrm = nfrm;

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");
    amx = dmatrx(1,ply1,1,ste1);
    amy = dmatrx(1,ply1,1,ste1);
    amz = dmatrx(1,ply1,1,ste1);
    bmx = dmatrx(1,ply2,1,ste2);
    bmy = dmatrx(1,ply2,1,ste2);
    bmz = dmatrx(1,ply2,1,ste2);
    xam = dmatrx(1,ply1,1,mon1);
    yam = dmatrx(1,ply1,1,mon1);
    zam = dmatrx(1,ply1,1,mon1);
    xbm = dmatrx(1,ply2,1,mon2);
    ybm = dmatrx(1,ply2,1,mon2);
    zbm = dmatrx(1,ply2,1,mon2);
    acx = dvectr(1,ply1);
    acy = dvectr(1,ply1);
    acz = dvectr(1,ply1);
    bcx = dvectr(1,ply2);
    bcy = dvectr(1,ply2);
    bcz = dvectr(1,ply2);

    //  PERFORM CALCULATIONS
    fprintf(stdout, "Calculations are being performed.\n");
    grga = 0.0;
    grba = 0.0;
    grea = 0.0;
    grgb = 0.0;
    grbb = 0.0;
    greb = 0.0;
    fram = fopen(pfrm_, "r");
    cord = fopen(pcrd_, "r");
    frb1 = fopen(prb1_, "w");
    frg1 = fopen(prg1_, "w");
    fre1 = fopen(pre1_, "w");
    frb2 = fopen(prb2_, "w");
    frg2 = fopen(prg2_, "w");
    fre2 = fopen(pre2_, "w");
    for(a=1; a<=nfrm; a++)
    {
        //  READ FRAME CONFIGURATION
        fprintf(stdout, "FRAME %5li of %5li\n", a, nfrm);
        fscanf(fram, "%li", &step);
        for(i=1; i<=ply1; i++)
        {
            for(j=1; j<=ste1; j++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                amx[i][j] = xdis;
                amy[i][j] = ydis;
                amz[i][j] = zdis;
            }
        }
        for(i=1; i<=ply2; i++)
        {
            for(j=1; j<=ste2; j++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                bmx[i][j] = xdis;
                bmy[i][j] = ydis;
                bmz[i][j] = zdis;
            }
        }

        //  COMPUTE MOLECULAR CENTER-OF-MASS COORDINATES
        for(i=1; i<=ply1; i++)
        {
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(j=1; j<=ste1; j++)
            {
                xdis = xdis + amx[i][j];
                ydis = ydis + amy[i][j];
                zdis = zdis + amz[i][j];
            }
            acx[i] = xdis/(double)ste1;
            acy[i] = ydis/(double)ste1;
            acz[i] = zdis/(double)ste1;
        }
        for(i=1; i<=ply2; i++)
        {
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(j=1; j<=ste2; j++)
            {
                xdis = xdis + bmx[i][j];
                ydis = ydis + bmy[i][j];
                zdis = zdis + bmz[i][j];
            }
            bcx[i] = xdis/(double)ste2;
            bcy[i] = ydis/(double)ste2;
            bcz[i] = zdis/(double)ste2;
        }

        //  COMPUTE MONOMER COORDINATES
        for(i=1; i<=ply1; i++)
        {
            for(j=1; j<=mon1; j++)
            {
                xdis = 0.0;
                ydis = 0.0;
                zdis = 0.0;
                for(k=1; k<=bds1; k++)
                {
                    xdis = xdis + amx[i][(j-1)*bds1+k];
                    ydis = ydis + amy[i][(j-1)*bds1+k];
                    zdis = zdis + amz[i][(j-1)*bds1+k];
                }
                xam[i][j] = xdis/(double)bds1;
                yam[i][j] = ydis/(double)bds1;
                zam[i][j] = zdis/(double)bds1;
            }
        }
        for(i=1; i<=ply2; i++)
        {
            for(j=1; j<=mon2; j++)
            {
                xdis = 0.0;
                ydis = 0.0;
                zdis = 0.0;
                for(k=1; k<=bds2; k++)
                {
                    xdis = xdis + bmx[i][(j-1)*bds2+k];
                    ydis = ydis + bmy[i][(j-1)*bds2+k];
                    zdis = zdis + bmz[i][(j-1)*bds2+k];
                }
                xbm[i][j] = xdis/(double)bds2;
                ybm[i][j] = ydis/(double)bds2;
                zbm[i][j] = zdis/(double)bds2;
            }
        }

        //  COMPUTE MONOMER BOND LENGTH
        dist = 0.0;
        for(i=1; i<=ply1; i++)
        {
            for(j=1; j<mon1; j++)
            {
                xdis = xam[i][j+1] - xam[i][j];
                ydis = yam[i][j+1] - yam[i][j];
                zdis = zam[i][j+1] - zam[i][j];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)((mon1-1)*ply1);
        grba = grba + dist;
        fprintf(frb1, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, grba/(double)a);
        dist = 0.0;
        for(i=1; i<=ply2; i++)
        {
            for(j=1; j<mon2; j++)
            {
                xdis = xbm[i][j+1] - xbm[i][j];
                ydis = ybm[i][j+1] - ybm[i][j];
                zdis = zbm[i][j+1] - zbm[i][j];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)((mon2-1)*ply2);
        grbb = grbb + dist;
        fprintf(frb2, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, grbb/(double)a);

        //  COMPUTE MOLECULAR RADIUS OF GYRATION
        dist = 0.0;
        for(i=1; i<=ply1; i++)
        {
            for(j=1; j<=ste1; j++)
            {
                xdis = amx[i][j] - acx[i];
                ydis = amy[i][j] - acy[i];
                zdis = amz[i][j] - acz[i];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)(ply1*ste1);
        grga = grga + dist;
        fprintf(frg1, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, grga/(double)a);
        dist = 0.0;
        for(i=1; i<=ply2; i++)
        {
            for(j=1; j<=ste2; j++)
            {
                xdis = bmx[i][j] - bcx[i];
                ydis = bmy[i][j] - bcy[i];
                zdis = bmz[i][j] - bcz[i];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)(ply2*ste2);
        grgb = grgb + dist;
        fprintf(frg2, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, grgb/(double)a);

        //  COMPUTE MOLECULAR END-TO-END DISTANCE
        dist = 0.0;
        for(i=1; i<=ply1; i++)
        {
            xdis = xam[i][1] - xam[i][mon1];
            ydis = yam[i][1] - yam[i][mon1];
            zdis = zam[i][1] - zam[i][mon1];
            dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
        }
        dist = dist/(double)(ply1);
        grea = grea + dist;
        fprintf(fre1, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, grea/(double)a);
        dist = 0.0;
        for(i=1; i<=ply2; i++)
        {
            xdis = xbm[i][1] - xbm[i][mon2];
            ydis = ybm[i][1] - ybm[i][mon2];
            zdis = zbm[i][1] - zbm[i][mon2];
            dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
        }
        dist = dist/(double)(ply2);
        greb = greb + dist;
        fprintf(fre2, "% 10li\t% 20.16E\t% 20.16E\n", step, dist, greb/(double)a);
    }
    fclose(fram);
    fclose(cord);
    fclose(frb1);
    fclose(frg1);
    fclose(fre1);
    fclose(frb2);
    fclose(frg2);
    fclose(fre2);

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
