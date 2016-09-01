//  ROUTINE CHECKS FOR TRAJECTORY JUMPS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutility.h"
#define MXSZ 500
void file_chck(FILE **, char *);
int main()
{
    char pfrm_[MXSZ], pcrd_[MXSZ];
    double *icmx, *icmy, *icmz, *ocmx, *ocmy, *ocmz,
        xdis, ydis, zdis, boxh, boxs;
    long i, j, m, t, nfrm, poly, mono;
    FILE *cord, *fram, *chck;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T448_66/";              //  INPUT: DATA DIRECTORY
    char cord_[] = "jcor_ETHY_T448_66";               //  INPUT: SITE COORDINATES
    char fram_[] = "fram_ETHY_T448_66";               //  INPUT: SIMULATION TIME STEPS
    char chck_[] = "jmps_ETHY_T448_66";               //  OUTPUT: BINARY SITE COORDINATES
    poly = 100;                                       //  NUMBER OF CHAINS
    mono = 66;                                        //  NUMBER OF SITES PER POLYMER
    boxs = 58.553;                                    //  MEAN BOX DIMENSION

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pfrm_, "%s%s", base_, fram_);

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
    icmx = dvectr(1,poly);
    icmy = dvectr(1,poly);
    icmz = dvectr(1,poly);
    ocmx = dvectr(1,poly);
    ocmy = dvectr(1,poly);
    ocmz = dvectr(1,poly);

    //  CHECK FOR TRAJECTORY JUMPS
    fprintf(stdout, "Trajectory jumps are being checked.\n");
    boxh = boxs/2.0;
    cord = fopen(pcrd_, "r");
    chck = fopen(chck_, "w");
    for(t=1; t<=nfrm; t++)
    {
        fprintf(stdout, "CYCLE: %5li of %5li\n", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            icmx[m] = 0.0;
            icmy[m] = 0.0;
            icmz[m] = 0.0;
            for(i=1; i<=mono; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &j, &xdis, &ydis, &zdis);
                icmx[m] = icmx[m] + xdis;
                icmy[m] = icmy[m] + ydis;
                icmz[m] = icmz[m] + zdis;
            }
            icmx[m] = icmx[m]/(double)mono;
            icmy[m] = icmy[m]/(double)mono;
            icmz[m] = icmz[m]/(double)mono;
            //  CHECK FOR IRREGULAR DISPLACEMENTS
            if(t>1)
            {
                xdis = fabs(ocmx[m] - icmx[m]);
                if(xdis > boxh)
                {
                    fprintf(chck, "X-VIOLATION: FRAME %5li\t CHAIN %5li\n", t, m);
                }
                ydis = fabs(ocmy[m] - icmy[m]);
                if(ydis > boxh)
                {
                    fprintf(chck, "Y-VIOLATION: FRAME %5li\t CHAIN %5li\n", t, m);
                }
                zdis = fabs(ocmz[m] - icmz[m]);
                if(zdis > boxh)
                {
                    fprintf(chck, "Z-VIOLATION: FRAME %5li\t CHAIN %5li\n", t, m);
                }
            }
            ocmx[m] = icmx[m];
            ocmy[m] = icmy[m];
            ocmz[m] = icmz[m];
        }
    }
    fclose(chck);
    fclose(cord);

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
