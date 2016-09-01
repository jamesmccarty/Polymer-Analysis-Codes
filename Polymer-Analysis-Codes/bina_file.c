//  ROUTINE CREATES A BINARY FILE OF MONOMER AND CENTER OF MASS COORDINATES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
    char pfrm_[MXSZ], pcrd_[MXSZ], pbmn_[MXSZ], pbcm_[MXSZ];
    double xdis, ydis, zdis, xadd, yadd, zadd;
    long i, j, m, t, nfrm, poly, mono, lsze;
    FILE *cord, *fram, *bmon, *bcom;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_HHPR_T453_96/";              //  INPUT: DATA DIRECTORY
    char cord_[] = "cord_HHPR_T453_96";               //  INPUT: SITE COORDINATES
    char fram_[] = "fram_HHPR_T453_96";               //  INPUT: SIMULATION TIME STEPS
    char bmon_[] = "bmon_HHPR_T453_96";               //  OUTPUT: BINARY SITE COORDINATES
    char bcom_[] = "bcom_HHPR_T453_96";               //  OUTPUT: BINARY COM COORDINATES
    poly = 1600;                                      //  NUMBER OF CHAINS
    mono = 96;                                        //  NUMBER OF SITES PER POLYMER

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pbmn_, "%s%s", base_, bmon_);
    sprintf(pbcm_, "%s%s", base_, bcom_);

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

    //  GENERATE BINARY COORDINATES FILE
    fprintf(stdout, "Binary coordinates are being generated.\n");
    lsze = sizeof(struct srcv);
    cord = fopen(pcrd_, "r");
    bmon = fopen(pbmn_, "w");
    bcom = fopen(pbcm_, "w");
    for(t=1; t<=nfrm; t++)
    {
        for(m=1; m<=poly; m++)
        {
            xadd = 0.0;
            yadd = 0.0;
            zadd = 0.0;
            for(i=1; i<=mono; i++)
            {
                fscanf(cord, "%li %lf %lf %lf", &j, &xdis, &ydis, &zdis);
                cont.x = xdis;
                cont.y = ydis;
                cont.z = zdis;
                fwrite(&cont, lsze, 1, bmon);
                xadd = xadd + xdis;
                yadd = yadd + ydis;
                zadd = zadd + zdis;
            }
            xadd = xadd/(double)mono;
            yadd = yadd/(double)mono;
            zadd = zadd/(double)mono;
            cont.x = xadd;
            cont.y = yadd;
            cont.z = zadd;
            fwrite(&cont, lsze, 1, bcom);
        }
    }
    fclose(bmon);
    fclose(bcom);
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
