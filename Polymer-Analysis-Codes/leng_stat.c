//  ROUTINE COMPUTES MOLECULAR LENGTH STATISTICS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutility.h"
#define MXSZ 500
void file_chck(FILE **, char *);
int main()
{
    char comd_[MXSZ], pcrd_[MXSZ], pfrm_[MXSZ], prbd_[MXSZ], prgc_[MXSZ], pret_[MXSZ];
    double **xste, **yste, **zste, **xmon, **ymon, **zmon,
        *xcom, *ycom, *zcom, xdis, ydis, zdis,
        dist, grgc, grgm, grbd, gret;
    long a, i, j, k, nfrm, poly, mono, bead, site, step;
    FILE *cord, *fram, *frbd, *frgc, *fret;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T450_0100/";            //  INPUT: DATA DIRECTORY
    char cord_[] = "cord_ETHY_T450_0100";             //  INPUT: SITE COORDINATES FILE
    char fram_[] = "fram_ETHY_T450_0100";             //  INPUT: SIMULATION TIME STEPS
    char dirn_[] = "LENG_ETHY_T450_0100/";            //  OUTPUT: RESULTS DIRECTORY
    char frbd_[] = "rbnd_ETHY_T450_0100";             //  OUTPUT: MONOMER BOND LENGTH
    char frgc_[] = "rgcm_ETHY_T450_0100";             //  OUTPUT: MOLECULAR GYRATION RADIUS
    char fret_[] = "rete_ETHY_T450_0100";             //  OUTPUT: POLYMER EXTENSION
    bead = 1;                                         //  NUMBER OF SITES IN MONOMER
    site = 100;                                       //  NUMBER OF SITES IN CHAIN
    poly = 48;                                        //  NUMBER OF CHAINS
    mono = site/bead;                                 //  NUMBER OF MONOMERS IN CHAIN

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(prbd_, "%s%s", dirn_, frbd_);
    sprintf(prgc_, "%s%s", dirn_, frgc_);
    sprintf(pret_, "%s%s", dirn_, fret_);
    sprintf(comd_, "mkdir %s", dirn_);
    system(comd_);
    frbd = fopen(prbd_, "w");
    frgc = fopen(prgc_, "w");
    fret = fopen(pret_, "w");

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
    xste = dmatrx(1,poly,1,site);
    yste = dmatrx(1,poly,1,site);
    zste = dmatrx(1,poly,1,site);
    xmon = dmatrx(1,poly,1,mono);
    ymon = dmatrx(1,poly,1,mono);
    zmon = dmatrx(1,poly,1,mono);
    xcom = dvectr(1,poly);
    ycom = dvectr(1,poly);
    zcom = dvectr(1,poly);

    //  PERFORM CALCULATIONS
    fprintf(stdout, "Calculations are being performed.\n");
    grbd = 0.0;
    grgm = 0.0;
    grgc = 0.0;
    gret = 0.0;
    fram = fopen(pfrm_, "r");
    cord = fopen(pcrd_, "r");
    for(a=1; a<=nfrm; a++)
    {
        //  READ FRAME CONFIGURATION
        fprintf(stdout, "FRAME %5li of %5li\n", a, nfrm);
        fscanf(fram, "%li", &step);
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=site; j++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                xste[i][j] = xdis;
                yste[i][j] = ydis;
                zste[i][j] = zdis;
            }
        }

        //  COMPUTE MOLECULAR CENTER-OF-MASS COORDINATES
        for(i=1; i<=poly; i++)
        {
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(j=1; j<=site; j++)
            {
                xdis = xdis + xste[i][j];
                ydis = ydis + yste[i][j];
                zdis = zdis + zste[i][j];
            }
            xcom[i] = xdis/(double)site;
            ycom[i] = ydis/(double)site;
            zcom[i] = zdis/(double)site;
        }

        //  COMPUTE MONOMER COORDINATES
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=mono; j++)
            {
                xdis = 0.0;
                ydis = 0.0;
                zdis = 0.0;
                for(k=1; k<=bead; k++)
                {
                    xdis = xdis + xste[i][(j-1)*bead+k];
                    ydis = ydis + yste[i][(j-1)*bead+k];
                    zdis = zdis + zste[i][(j-1)*bead+k];
                }
                xmon[i][j] = xdis/(double)bead;
                ymon[i][j] = ydis/(double)bead;
                zmon[i][j] = zdis/(double)bead;
            }
        }

        //  COMPUTE MONOMER BOND LENGTH
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<mono; j++)
            {
                xdis = xmon[i][j+1] - xmon[i][j];
                ydis = ymon[i][j+1] - ymon[i][j];
                zdis = zmon[i][j+1] - zmon[i][j];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)((mono-1)*poly);
        grbd = grbd + dist;
        fprintf(frbd, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, grbd/(double)a);

        //  COMPUTE MOLECULAR RADIUS OF GYRATION
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=site; j++)
            {
                xdis = xste[i][j] - xcom[i];
                ydis = yste[i][j] - ycom[i];
                zdis = zste[i][j] - zcom[i];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)(poly*site);
        grgc = grgc + dist;
        fprintf(frgc, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, grgc/(double)a);

        //  COMPUTE MOLECULAR END-TO-END DISTANCE
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            xdis = xmon[i][1] - xmon[i][mono];
            ydis = ymon[i][1] - ymon[i][mono];
            zdis = zmon[i][1] - zmon[i][mono];
            dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
        }
        dist = dist/(double)(poly);
        gret = gret + dist;
        fprintf(fret, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, gret/(double)a);
    }
    fclose(fram);
    fclose(cord);
    fclose(frbd);
    fclose(frgc);
    fclose(fret);

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
