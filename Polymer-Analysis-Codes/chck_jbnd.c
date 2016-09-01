//  ROUTINE COMPUTES REJOINED SITE COORDINATES USING PERIODIC BOUNDARY CONDITIONS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MXTH 2.0
#define MNTH 1.0
#define MXSZ 500
void file_chck(FILE **, char *);
int main()
{
    char pfrm_[MXSZ], pcrd_[MXSZ];
    double boxs, boxh, newx, newy, newz, oewx, oewy, oewz;
    long i, m, t, id, poly, mono, nfrm;
    FILE *cord, *fram, *chck;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T448_66/";              //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T448_66";               //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "cord_ETHY_T448_66";               //  INPUT: SITE COORDINATES
    char chck_[] = "jbnd_ETHY_T448_66";               //  OUTPUT: ENSEMBLE COORDINATES
    mono = 66;                                        //  NUMBER OF SITES PER POLYMER
    poly = 100;                                       //  NUMBER OF CHAINS
    boxs = 58.553;

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

    //  JOIN BROKEN BONDS
    fprintf(stdout, "Broken bonds are being checked and joined.\n");
    cord = fopen(pcrd_, "r");
    chck = fopen(chck_, "w");
    boxh = boxs/2.0;
    for(t=1; t<=nfrm; t++)
    {
        fprintf(stdout, "CYCLE: % 5li of % 5li\n", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            fscanf(cord, "%li %lf %lf %lf", &id, &oewx, &oewy, &oewz);
            fprintf(chck, "% 5li\t% 12.6lf\t% 12.6lf\t% 12.6lf\n", id, oewx, oewy, oewz);
            for(i=1; i<=mono-1; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &id, &newx, &newy, &newz);
                while(fabs(oewx - newx)>=boxh)
                {
                    if((oewx-newx)<0)
                    {
                        newx = newx - boxs;
                    }
                    else
                    {
                        newx = newx + boxs;
                    }
                }
                while(fabs(oewy - newy)>=boxh)
                {
                    if((oewy-newy)<0)
                    {
                        newy = newy - boxs;
                    }
                    else
                    {
                        newy = newy + boxs;
                    }
                }
                while(fabs(oewz - newz)>=boxh)
                {
                    if((oewz-newz)<0)
                    {
                        newz = newz - boxs;
                    }
                    else
                    {
                        newz = newz + boxs;
                    }
                }
                oewx = newx;
                oewy = newy;
                oewz = newz;
                fprintf(chck, "% 5li\t% 12.6lf\t% 12.6lf\t% 12.6lf\n", id, oewx, oewy, oewz);
            }
        }
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
