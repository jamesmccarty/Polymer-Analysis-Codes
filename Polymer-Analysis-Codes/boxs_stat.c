//  ROUTINE COMPUTES BOX SIZE STATISTICS (FOR NPT ENSEMBLES)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define BXMN 0.0
#define BXMX 1000.0
#define MXSZ 500
void file_chck(FILE **, char *);
int main()
{
    char comd_[MXSZ], pfrm_[MXSZ], pbox_[MXSZ], pbav_[MXSZ], psta_[MXSZ];
    double xsde, ysde, zsde, mxsz, mnsz, mean, sdev;
    long i, nfrm;
    FILE *boxs, *fram, *fbav, *fsta;

    //  INITIALIZE ALGORITHM PARAMETERS
    char base_[] = "DATA_ETHY_T450_0100/";            //  INPUT: DATA DIRECTORY
    char boxs_[] = "boxs_ETHY_T450_0100";             //  INPUT: BOX DIMENSIONS FILE
    char fram_[] = "fram_ETHY_T450_0100";             //  INPUT: FRAME STAMP FILE
    char dirn_[] = "LENG_ETHY_T450_0100/";            //  OUTPUT: RESULTS DIRECTORY
    char fbav_[] = "cbav_ETHY_T450_0100";             //  OUTPUT: BOX SIZE AVERAGES
    char fsta_[] = "bsta_ETHY_T450_0100";             //  OUTPUT: BOX SIZE STATISTICS

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pbox_, "%s%s", base_, boxs_);
    sprintf(pbav_, "%s%s", dirn_, fbav_);
    sprintf(psta_, "%s%s", dirn_, fsta_);
    sprintf(comd_, "mkdir %s", dirn_);
    system(comd_);

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&boxs, pbox_);
    fprintf(stdout, "Status of input files is satisfactory.\n");

    //  COUNT SIMULATION FRAMES
    fprintf(stdout, "Frames are being counted.\n");
    fram = fopen(pfrm_, "r");
    nfrm = 0;
    while(fscanf(fram, "%li", &i)!=EOF)
    {
        nfrm++;
    }
    fclose(fram);
    nfrm = nfrm;

    //  COMPUTE BOX MEAN DIMENSION
    fprintf(stdout, "Box mean dimension is being computed.\n");
    mxsz = BXMN;
    mnsz = BXMX;
    mean = 0.0;
    boxs = fopen(pbox_, "r");
    fbav = fopen(pbav_, "w");
    fsta = fopen(psta_, "w");
    for(i=1; i<=nfrm; i++)
    {
        fscanf(boxs, "%lf %lf %lf", &xsde, &ysde, &zsde);
        mean = mean + xsde;
        if(xsde < mnsz)
        {
            mnsz = xsde;
        }
        if(xsde > mxsz)
        {
            mxsz = xsde;
        }
        fprintf(fbav, "% 20.16E\t% 20.16E\n",
            xsde, mean/(double)i);
    }
    mean = mean/(double)(nfrm);
    fclose(fbav);

    //  COMPUTE STANDARD DEVIATION
    fprintf(stdout, "Box dimension standard deviation is being computed.\n");
    sdev = 0.0;
    rewind(boxs);
    while(fscanf(boxs, "%lf %lf %lf", &xsde, &ysde, &zsde)!=EOF)
    {
        sdev = sdev + (xsde-mean)*(xsde-mean);
    }
    fclose(boxs);
    sdev = sqrt(sdev/(double)(nfrm));

    //  PRINT BOX STATISTICS SUMMARY
    fprintf(stdout, "Box statistics summary is being printed.\n");
    fprintf(fsta, "MAXIMUM: % 20.16E\n", mxsz);
    fprintf(fsta, "MINIMUM: % 20.16E\n", mnsz);
    fprintf(fsta, "AVERAGE: % 20.16E\n", mean);
    fprintf(fsta, "STD DEV: % 20.16E\n", sdev);
    fclose(fsta);

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
