//  ROUTINE COMPUTES MEAN SQUARE DISPLACEMENT AND SECOND-ORDER NON-GAUSSIAN
//  ALPHA PARAMETER FOR MONOMER AND CENTER OF MASS COORDINATES
#include <stdio.h>
#include <stdlib.h>
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
    char comd_[MXSZ], pcrd_[MXSZ], pfrm_[MXSZ], pmsm_[MXSZ], pmsc_[MXSZ];
    double xdis, ydis, zdis, xacu, yacu, zacu, alfa, prod, fact,
        *xcor, *ycor, *zcor, *dsqr, *dqrt;
    long i, m, t1, t2, nfrm, poly, mono, *step, lsze, fsze;
    FILE *cord, *fram, *msdm, *msdc;

    //  INITIALIZE ALGORITHM VARIABLES
    char base_[] = "DATA_ETHY_T450_0100/";            //  INPUT: DATA DIRECTORY
    char fram_[] = "fram_ETHY_T450_0100";             //  INPUT: SIMULATION TIME STEPS
    char cord_[] = "bmon_ETHY_T450_0100";             //  INPUT: MONOMER COORDINATES IN BINARY
    char dirn_[] = "MDIS_ETHY_T450_0100/";            //  OUTPUT: RESULTS DIRECTORY
    char msdm_[] = "msdm_ETHY_T450_0100";             //  OUTPUT: MONOMER DISPLACEMENT/ALPHA
    char msdc_[] = "msdc_ETHY_T450_0100";             //  OUTPUT: COM DISPLACEMENT/ALPHA
    poly = 48;                                        //  NUMBER OF POLYMERS
    mono = 100;                                       //  NUMBER OF MONOMERS PER POLYMER
    fact = 0.000001;                                  //  TIME CONVERSION FACTOR

    //  OPEN PROGRAM SEQUENCE
    system("clear");
    setlinebuf(stdout);
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pmsm_, "%s%s", dirn_, msdm_);
    sprintf(pmsc_, "%s%s", dirn_, msdc_);
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
    fprintf(stdout, "System resources are being allocated.\n");
    step = lvectr(1, nfrm);
    xcor = dvectr(1, nfrm);
    ycor = dvectr(1, nfrm);
    zcor = dvectr(1, nfrm);
    dsqr = dvectr(1, nfrm);
    dqrt = dvectr(1, nfrm);
    for(i=1; i<=nfrm; i++)
    {
        dsqr[i] = 0.0;
        dqrt[i] = 0.0;
    }
    lsze = sizeof(struct srcv);
    fsze = mono*poly*lsze;

    //  READ STEPS FILE
    fprintf(stdout, "Frames are being read.\n");
    rewind(fram);
    for(i=1; i<=nfrm; i++)
    {
        fscanf(fram, "%li", &step[i]);
    }
    fclose(fram);

    //  CALCULATE MONOMER DISPLACEMENT
    fprintf(stdout, "Monomer displacement is being calculated.\n");
    cord = fopen(pcrd_, "r");
    for(m=1; m<=poly; m++)
    {
        fprintf(stdout, "MONOMER CYCLE: %5li of %5li\t", m, poly);
        si = time(NULL);
        for(i=1; i<=mono; i++)
        {
            for(t1=1; t1<=nfrm; t1++)
            {
                if(t1==1)
                {
                    fseek(cord, ((m-1)*mono+(i-1))*lsze, SEEK_SET);
                }
                fread(&cont, lsze, 1, cord);
                xcor[t1] = cont.x;
                ycor[t1] = cont.y;
                zcor[t1] = cont.z;
                fseek(cord, fsze-lsze, SEEK_CUR);
            }
            for(t1=1; t1<=nfrm; t1++)
            {
                for(t2=t1; t2<=nfrm; t2++)
                {
                    xdis = xcor[t2] - xcor[t1];
                    ydis = ycor[t2] - ycor[t1];
                    zdis = zcor[t2] - zcor[t1];
                    prod = xdis*xdis + ydis*ydis + zdis*zdis;
                    dsqr[t2-t1+1] = dsqr[t2-t1+1] + prod;
                    dqrt[t2-t1+1] = dqrt[t2-t1+1] + prod*prod;
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }

    //  PRINT RESULTS
    fprintf(stdout, "Results are being printed.\n");
    msdm = fopen(pmsm_, "w");
    for(t1=2; t1<=nfrm; t1++)
    {
        dsqr[t1] = dsqr[t1]/(double)(poly*mono*(nfrm-t1+1));
        dqrt[t1] = dqrt[t1]/(double)(poly*mono*(nfrm-t1+1));
        alfa = 3.0/5.0*dqrt[t1]/(dsqr[t1]*dsqr[t1]) - 1.0;
        fprintf(msdm, "% 20.16E\t% 20.16E\t% 20.16E\t% 20.16E\n",
            fact*(double)step[t1], dsqr[t1], dqrt[t1], alfa);
    }
    fclose(msdm);

    //  CALCULATE CENTER OF MASS DISPLACEMENT
    for(i=1; i<=nfrm; i++)
    {
        dsqr[i] = 0.0;
        dqrt[i] = 0.0;
    }
    fprintf(stdout, "Center of mass displacement is being calculated.\n");
    rewind(cord);
    for(m=1; m<=poly; m++)
    {
        fprintf(stdout, "CENTROID CYCLE: %5li of %5li\t", m, poly);
        si = time(NULL);
        for(t1=1; t1<=nfrm; t1++)
        {
            if(t1==1)
            {
                fseek(cord, ((m-1)*mono)*lsze, SEEK_SET);
            }
            xacu = 0.0;
            yacu = 0.0;
            zacu = 0.0;
            for(i=1; i<=mono; i++)
            {
                fread(&cont, lsze, 1, cord);
                xacu = xacu + cont.x;
                yacu = yacu + cont.y;
                zacu = zacu + cont.z;
            }
            fseek(cord, fsze-lsze*mono, SEEK_CUR);
            xcor[t1] = xacu/(double)mono;
            ycor[t1] = yacu/(double)mono;
            zcor[t1] = zacu/(double)mono;
        }
        for(t1=1; t1<=nfrm; t1++)
        {
            for(t2=t1; t2<=nfrm; t2++)
            {
                xdis = xcor[t2] - xcor[t1];
                ydis = ycor[t2] - ycor[t1];
                zdis = zcor[t2] - zcor[t1];
                prod = xdis*xdis + ydis*ydis + zdis*zdis;
                dsqr[t2-t1+1] = dsqr[t2-t1+1] + prod;
                dqrt[t2-t1+1] = dqrt[t2-t1+1] + prod*prod;
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  PRINT RESULTS
    fprintf(stdout, "Results are being printed.\n");
    msdc = fopen(pmsc_, "w");
    for(t1=2; t1<=nfrm; t1++)
    {
        dsqr[t1] = dsqr[t1]/(double)(poly*(nfrm-t1+1));
        dqrt[t1] = dqrt[t1]/(double)(poly*(nfrm-t1+1));
        alfa = 3.0/5.0*dqrt[t1]/(dsqr[t1]*dsqr[t1]) - 1.0;
        fprintf(msdc, "% 20.16E\t% 20.16E\t% 20.16E\t% 20.16E\n",
            step[t1]*fact, dsqr[t1], dqrt[t1], alfa);
    }
    fclose(msdc);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has terminated successfully.\n");
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
