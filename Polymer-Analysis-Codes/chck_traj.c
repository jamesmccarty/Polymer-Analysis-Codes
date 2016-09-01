//  ROUTINE REGISTERS TRAJECTORIES OF PARTICLES IN ENSEMBLE
#include <stdio.h>
#include <stdlib.h>
#include "nrutility.h"
#define MXSZ 500
void file_chck(FILE **, char *);
struct srcv
{
    double x, y, z;
};
int main()
{
    struct srcv cont;
    char pfrm_[MXSZ], pcrd_[MXSZ], pchk_[MXSZ], pvar_[MXSZ];
    double *xcor, *ycor, *zcor;
    long i, m, t, poly, mono, nfrm, lsze, fsze;
    FILE *cord, *fram, *chck;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T448_66/";            //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T448_66";             //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "jtrj_ETHY_T448_66";             //  INPUT: SITE COORDINATES IN BINARY
    char chck_[] = "traj_ETHY_T448_66";             //  OUTPUT: TRAJECTORY FILE
    mono = 66;                                        //  NUMBER OF SITES PER POLYMER
    poly = 100;                                       //  NUMBER OF CHAINS

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);

    //  QUERY STATS OF INPUT FILES
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
    xcor = dvectr(1,nfrm);
    ycor = dvectr(1,nfrm);
    zcor = dvectr(1,nfrm);
    lsze = sizeof(struct srcv);
    fsze = mono*poly*lsze;

    //  COLLECT PARTICLE TRAJECTORIES
    fprintf(stdout, "Particle trajectories are being checked.\n");
    cord = fopen(pcrd_, "r");
    for(m=1; m<=poly; m++)
    {
        fprintf(stdout, "CYCLE: % 5li of % 5li\n", m, poly);
        for(t=1; t<=nfrm; t++)
        {
            xcor[t] = 0.0;
            ycor[t] = 0.0;
            zcor[t] = 0.0;
        }
        for(i=1; i<=mono; i++)
        {
            sprintf(pchk_, "%s_%li_%li", chck_, m, i);
            chck = fopen(pchk_, "w");
            for(t=1; t<=nfrm; t++)
            {
                if(t==1)
                {
                    fseek(cord, ((m-1)*mono+(i-1))*lsze, SEEK_SET);
                }
                fread(&cont, lsze, 1, cord);
                xcor[t] = xcor[t] + cont.x;
                ycor[t] = ycor[t] + cont.y;
                zcor[t] = zcor[t] + cont.z;
                fseek(cord, fsze-lsze, SEEK_CUR);
                fprintf(chck, "% 5li\t% 20.10lf\t% 20.10lf\t% 20.10lf\n",
                t, cont.x, cont.y, cont.z);
            }
            fclose(chck);
        }
        sprintf(pvar_, "%s_%li_%li", chck_, m, (long)0);
        chck = fopen(pvar_, "w");
        for(t=1; t<=nfrm; t++)
        {
            xcor[t] = xcor[t]/(double)(mono);
            ycor[t] = ycor[t]/(double)(mono);
            zcor[t] = zcor[t]/(double)(mono);
            fprintf(chck, "% 5li\t% 20.10lf\t% 20.10lf\t% 20.10lf\n",
            t, xcor[t], ycor[t], zcor[t]);
        }
        fclose(chck);
    }
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
