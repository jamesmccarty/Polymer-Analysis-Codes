//  ROUTINE REGISTERS TRAJECTORIES OF PARTICLES IN ENSEMBLE
#include <stdio.h>
#include <stdlib.h>
#include "nrutility.h"
#define MXSZ 500
long nrst_long(double);
void file_chck(FILE **, char *);
struct srcv
{
    double x, y, z;
};
int main()
{
    struct srcv cont;
    char pfrm_[MXSZ], pcrd_[MXSZ];
    double **xcor, **ycor, **zcor, **xraw, **yraw, **zraw,
        dx, dy, dz, boxs;
    long i, m, t, u, poly, mono, nfrm, lsze, fsze;
    FILE *cord, *fram, *chck, *ccor;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T448_66/";              //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T448_66";               //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "bmon_ETHY_T448_66";               //  INPUT: SITE COORDINATES IN BINARY
    char chck_[] = "jtrj_ETHY_T448_66";               //  OUTPUT: TRAJECTORY FILE
    char ccor_[] = "jcor_ETHY_T448_66";               //  OUTPUT: TRAJECTORY FILE
    mono = 66;                                        //  NUMBER OF SITES PER POLYMER
    poly = 100;                                       //  NUMBER OF CHAINS
    boxs = 58.553;                                    //  LINEAR DIMENSION OF SIMULATION BOX

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
    xcor = dmatrx(1,nfrm,1,poly);
    ycor = dmatrx(1,nfrm,1,poly);
    zcor = dmatrx(1,nfrm,1,poly);
    xraw = dmatrx(1,nfrm,1,poly);
    yraw = dmatrx(1,nfrm,1,poly);
    zraw = dmatrx(1,nfrm,1,poly);
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
            xraw[t][m] = 0.0;
            yraw[t][m] = 0.0;
            zraw[t][m] = 0.0;
        }
        for(i=1; i<=mono; i++)
        {
            for(t=1; t<=nfrm; t++)
            {
                if(t==1)
                {
                    fseek(cord, ((m-1)*mono+(i-1))*lsze, SEEK_SET);
                }
                fread(&cont, lsze, 1, cord);
                xraw[t][m] = xraw[t][m] + cont.x;
                yraw[t][m] = yraw[t][m] + cont.y;
                zraw[t][m] = zraw[t][m] + cont.z;
                fseek(cord, fsze-lsze, SEEK_CUR);
            }
        }
        for(t=1; t<=nfrm; t++)
        {
            xraw[t][m] = xraw[t][m]/(double)(mono);
            yraw[t][m] = yraw[t][m]/(double)(mono);
            zraw[t][m] = zraw[t][m]/(double)(mono);
            xcor[t][m] = xraw[t][m];
            ycor[t][m] = yraw[t][m];
            zcor[t][m] = zraw[t][m];
        }
        for(t=1; t<nfrm; t++)
        {
            dx = xraw[t+1][m] - xraw[t][m];
            dy = yraw[t+1][m] - yraw[t][m];
            dz = zraw[t+1][m] - zraw[t][m];
            for(u=t+1; u<=nfrm; u++)
            {
                xcor[u][m] = xcor[u][m] - boxs*(double)nrst_long(dx/boxs);
                ycor[u][m] = ycor[u][m] - boxs*(double)nrst_long(dy/boxs);
                zcor[u][m] = zcor[u][m] - boxs*(double)nrst_long(dz/boxs);
            }
        }
    }
    fclose(cord);

    //  APPLY CORRECTIONS TO COORDINATES
    fprintf(stdout, "Particle trajectories are being corrected.\n");
    cord = fopen(pcrd_, "r");
    chck = fopen(chck_, "w");
    ccor = fopen(ccor_, "w");
    for(t=1; t<=nfrm; t++)
    {
        fprintf(stdout, "CYCLE: % 5li of % 5li\n", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
                fread(&cont, lsze, 1, cord);
                cont.x = cont.x - xraw[t][m] + xcor[t][m];
                cont.y = cont.y - yraw[t][m] + ycor[t][m];
                cont.z = cont.z - zraw[t][m] + zcor[t][m];
                fwrite(&cont, lsze, 1, chck);
                fprintf(ccor, "% 5li\t% 12.6lf\t% 12.6lf\t% 12.6lf\n",
                    (long)1, cont.x, cont.y, cont.z);
            }
        }
    }
    fclose(cord);
    fclose(ccor);
    fclose(chck);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
    return(0);
}

//  NEAREST INTEGER HANDLER
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
