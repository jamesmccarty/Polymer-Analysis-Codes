//  ROUTINE COMPUTES INTRAMOLECULAR MONOMER/MONOMER
//  CORRELATION FUNCTION
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
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], pwok_[MXSZ];
    double **xcor, **ycor, **zcor, *wkit,
        kmax, kval, delk, xdis, ydis, zdis, dist;
    long i, j, k, m, t, knum, mono, poly, nfrm, lsze;
    FILE *fram, *cord, *wofk;

    //  INITIALIZE PARAMETERS HERE:
    char base_[] = "DATA_SNPR_T453_96/";              //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_SNPR_T453_96";               //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "bmon_SNPR_T453_96";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "HOFR_SNPR_T453_96/";              //  OUTPUT: RESULTS DIRECTORY
    char wofk_[] = "wmmk_SNPR_T453_96";               //  OUTPUT: INTRAMOLECULAR CORRELATION FUNCTION
    poly = 1600;                                      //  NUMBER OF CHAINS
    mono = 96;                                        //  NUMBER OF MONOMERS PER CHAIN
    delk = 0.01;                                      //  INCREMENTAL FACTOR OF WAVE VECTOR
    kmax = 2.0;                                      //  MAXIMUM MAGNITUDE OF WAVE VECTOR

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pwok_, "%s%s", dirn_, wofk_);
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
    while((fscanf(fram, "%li", &i))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);
    nfrm = nfrm;

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");
    knum = (long)(kmax/delk);
    xcor = dmatrx(1,poly,1,mono);
    ycor = dmatrx(1,poly,1,mono);
    zcor = dmatrx(1,poly,1,mono);
    wkit = dvectr(0,knum);
    for(k=0; k<=knum; k++)
    {
        wkit[k] = 0.0;
    }

    //  CALCULATE FORM FACTOR
    fprintf(stdout, "Form factor is being calculated.\n");
    lsze = sizeof(struct srcv);
    cord = fopen(pcrd_, "r");
    for(t=1; t<=nfrm; t++)
    {
        si = time(NULL);
        fprintf(stdout, "CYCLE %5li OF %5li\t", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
                fread(&cont, lsze, 1, cord);
                xcor[m][i] = cont.x;
                ycor[m][i] = cont.y;
                zcor[m][i] = cont.z;
            }
        }
        for(k=0; k<=knum; k++)
        {
            kval = (double)(k*delk);
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<mono; i++)
                {
                    for(j=i+1; j<=mono; j++)
                    {
                        xdis = xcor[m][i] - xcor[m][j];
                        ydis = ycor[m][i] - ycor[m][j];
                        zdis = zcor[m][i] - zcor[m][j];
                        dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                        if((kval*dist)!=0.0)
                        {
                            wkit[k] = wkit[k] + 2.0*sin(kval*dist)/(kval*dist);
                        }
                        else
                        {
                            wkit[k] = wkit[k] + 2.0;
                        }
                    }
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR FORM FACTOR
    fprintf(stdout, "Results are being printed.\n");
    wofk = fopen(pwok_, "w");
    for(k=0; k<=knum; k++)
    {
        kval = (double)(k*delk);
        wkit[k] = wkit[k] + (double)(nfrm*poly*mono);
        fprintf(wofk, "% 20.16E\t% 20.16E\n", kval, wkit[k]/(double)(nfrm*poly*mono));
    }
    fclose(wofk);

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
