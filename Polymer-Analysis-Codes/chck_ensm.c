//  ROUTINE COMPUTES CENTER OF MASS COORDINATE OF ENSEMBLE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MXSZ
void file_chck(FILE **, char *);
int main()
{
    char pfrm_[MXSZ], pcrd_[MXSZ];
    double corx, cory, corz, ensx, ensy, ensz, boxs, boxh;
    long i, m, t, id, poly, mono, nfrm;
    FILE *cord, *fram, *chck;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T448_66/";              //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T448_66";               //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "cord_ETHY_T448_66";               //  INPUT: SITE COORDINATES
    char chck_[] = "ensm_ETHY_T448_66";               //  OUTPUT: ENSEMBLE COORDINATES
    mono = 66;                                        //  NUMBER OF SITES PER POLYMER
    poly = 100;                                       //  NUMBER OF CHAINS
    boxs = 58.553;                                    //  LINEAR DIMENSIONS OF SIMULATION BOX

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

    //  CALCULATION CENTER OF MASS COORDINATE OF ENSEMBLE
    fprintf(stdout, "Ensemble center of mass coordinate is being calculated.\n");
    cord = fopen(pcrd_, "r");
    chck = fopen(chck_, "w");
    boxh = boxs/2.0;
    for(t=1; t<=nfrm; t++)
    {
        fprintf(stdout, "CYCLE: % 5li of % 5li\n", t, nfrm);
        ensx = 0.0;
        ensy = 0.0;
        ensz = 0.0;
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &id, &corx, &cory, &corz);
                if(corx>=boxh)
                {
                    corx = corx - boxs*(double)((long)(corx/boxs));
                }
                if(corx<=-boxh)
                {
                    corx = corx + boxs*(double)((long)(corx/boxs));
                }
                if(cory>=boxh)
                {
                    cory = cory - boxs*(double)((long)(cory/boxs));
                }
                if(cory<=-boxh)
                {
                    cory = cory + boxs*(double)((long)(cory/boxs));
                }
                if(corz>=boxh)
                {
                    corz = corz - boxs*(double)((long)(corz/boxs));
                }
                if(corz<=-boxh)
                {
                    corz = corz + boxs*(double)((long)(corz/boxs));
                }

                ensx = ensx + corx;
                ensy = ensy + cory;
                ensz = ensz + corz;
            }
        }
        ensx = ensx/(double)(poly*mono);
        ensy = ensy/(double)(poly*mono);
        ensz = ensz/(double)(poly*mono);
        fprintf(chck, "% 5li\t% 20.16lf\t% 20.16lf\t% 20.16lf\n",
            t, ensx, ensy, ensz);
    }
    fclose(cord);
    fclose(chck);

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
