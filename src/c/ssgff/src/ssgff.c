/*************************************************************************
*                                                                        *
*   Module: SSgff                                                        *
*                                                                        *
*   Formatted output of geneid (GFF, default and extended)               *
*                                                                        *
*   This file is part of the SSgff  Distribution                         *
*                                                                        *
*     Copyright (C) 2001 - Genis PARRA FARRE                             *
*                          Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           *
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

#include "ssgff.h"

/* SSgff setup flags */
int
/* sites to print */
    SFP = 0, SDP = 0, SAP = 0, STP = 0, ALLS = 0,
/* exons to print */
    EFP = 0, EIP = 0, ETP = 0, ESP = 0, ALLE = 0, INT = 0,
/* transcrip */
    CDS = 0, TRN = 0,
/* Verbose flag (memory/processing information) */
    VRB = 0;

/* crash debug dk 20150917
struct rusage r_usage;
*/

int main(int argc, char *argv[]) {

    /* input files names */
    char fasta_fn[FILENAMELENGTH], gff_fn[FILENAMELENGTH];

    /* sequence */
    char Locus[MAXLINE];
    char *Sequence;

    /* Auxiliary arrays */
    char *AuxExon;
    char *AuxGene;
    char *AuxExon_tmp;

    /* exons array */
    exonGFF   *exons;
    pack_gene *genes;

    /* variables */
    char out_string[MAXSTRING];
    long numex;
    long numgen;
    long long_cc, long_gg;
    long length, lengthCDS;
    long pos1, pos2;
    long pos1_site, pos2_site, pos1_intron, pos2_intron;

    /* Auxiliary pointers */
    exonGFF *AuxExonPtr;
    exonGFF *AuxExonPtr_tmp;

    /* dictionary variable for genes names */
    dict *dd_dict;

    /* input files */
    FILE *fasta_fh;
    FILE *gff_fh;

    /* crash debug
    getrusage(RUSAGE_SELF,&r_usage);
    printf(" 00 Memory usage = %ld\n",r_usage.ru_maxrss);
    crash debug end */

    /* memory for the sequence */
    if ((Sequence = (char *) calloc(SEQLENGTH, sizeof (char))) == NULL) {
        printf("Not enough memory: Sequence");
    }

    /* crash debug
    getrusage(RUSAGE_SELF,&r_usage);
    printf(" AAA 01 Memory usage = %ld\n",r_usage.ru_maxrss);
    crash debug end */

    if ((AuxExon
             = (char *) calloc(MAXAA * LENGTHCODON, sizeof (char))) == NULL) {
        printf("Not enough memory: AuxExon");
    }

    /* crash debug
    getrusage(RUSAGE_SELF,&r_usage);
    printf(" AAA 02 Memory usage = %ld\n",r_usage.ru_maxrss);
    crash debug end */

    if ((AuxGene
             = (char *) calloc(MAXAA * LENGTHCODON, sizeof (char))) == NULL) {
        printf("Not enough memory: AuxGene");
    }

    /* crash debug
    getrusage(RUSAGE_SELF,&r_usage);
    printf(" AAA 03 Memory usage = %ld\n",r_usage.ru_maxrss);
    crash debug end */

    if ((AuxExon_tmp
             = (char *) calloc(MAXAA * LENGTHCODON, sizeof (char))) == NULL) {
        printf("Not enough memory: AuxExon_tmp");
    }

    if ((exons = (exonGFF *) calloc(NUMEXONS, sizeof (exonGFF))) == NULL) {
        printf("Not enough memory: exons");
    }

    if ((genes = (pack_gene *) calloc(NUMEXONS, sizeof (pack_gene))) == NULL) {
        printf("Not enough memory: genes");
    }

    if ((AuxExonPtr_tmp = (exonGFF *) malloc(sizeof (exonGFF))) == NULL) {
        printf("Not enough memory: AuxExonPtr_tmp");
    }

    if ((AuxExonPtr = (exonGFF *) malloc(sizeof (exonGFF))) == NULL) {
        printf("Not enough memory: AuxExonPtr");
    }

    if ((dd_dict = (dict *) malloc(sizeof (dict))) == NULL) {
        printf("Not enough memory: dictionary of exon types");
    }

    /* crash debug
     getrusage(RUSAGE_SELF,&r_usage);
     printf(" XXX 01 Memory usage = %ld\n",r_usage.ru_maxrss);
     crash debug end */

    /*  Read setup options */
    readargv(argc, argv, fasta_fn, gff_fn);
    /* debug
    printf("\n\n\t\t\t** Executing ssgff 2015 **\n\n");
    */

    /* Open files */

    /* debug
    printf("fasta_fn: %s\n",  fasta_fn);
    printf("gff_fn:   %s\n",    gff_fn);
    */

    fasta_fh = OpenFile(fasta_fn);
    gff_fh   = OpenFile(gff_fn);

    /* read  the sequence */
    /* sprintf (out_string, "Reading FASTA file: %s", fasta_fn);
     printf("001 fasta_fn: %s\n",  fasta_fn);
    printf(out_string); #not working */
    /* printMess (out_string); */

    length = ReadSequence(fasta_fh, Locus, Sequence);
    /* debug
    printf("seq_len_fasta Locus: %lu %s \n",  length, Locus);
    */

    /*
    sprintf (out_string, "DNA sequence %s  readed (%ld bp)", Locus, length);
    printMess (out_string);
    */

    /* read exons */
    /*
    sprintf (out_string, "Reading GFF file: %s ", gff_fn);
    printMess (out_string);
    */
    /* debugg valgrind */
    numgen = 0;
    /* resetDict(dd_dict); */

    numex = ReadExonsGFF(gff_fh, exons, genes, dd_dict, &numgen, length);
    /*
    sprintf (out_string, "%ld exons readed grouped on %ld genes ", numex, numgen);
    printMess (out_string);
    */
    /* assert (numex   >=1);
       assert (numgen  >=1); */
    /* debug
    printf("num_ex: %lu\n",  numex);
    printf("num_gen: %lu\n",  numgen);
    */

    /* Scanning and printting the subsequences */
    for (long_gg = 0; long_gg < numgen; long_gg++) {

        /**** dumps core 20150917a */

        sprintf(out_string, "Extracting information from group:  %s , %ld exons",
                (genes + long_gg)->First->Group, (genes + long_gg)->numExons + 1);

        printMess(out_string);
        /**** dumps core 20150917a ****/

        lengthCDS = 0;

        /* Initializing firts exon before the loop */
        AuxExonPtr = (genes + long_gg)->First;

        for (long_cc = 0; long_cc <= (genes + long_gg)->numExons; long_cc++) {

            /* For other exons except the first one */
            if (long_cc != 0) {
                AuxExonPtr_tmp = AuxExonPtr;
                AuxExonPtr     = AuxExonPtr_tmp->NextExon;

                /* Intron setting part2 and Printing */
                if (INT) {
                    pos2_intron = (AuxExonPtr->Position1) - 1;
                    printHeadIntron(AuxExonPtr->Group, long_cc, Locus);
                    printFasta(pos1_intron, pos2_intron, AuxExonPtr->Strand,
                               Sequence);
                }
            }

            /* Setting intron pos1 for the next loop */
            pos1_intron = (AuxExonPtr->Position2) + 1;

            /* Defining exon position1 and position2 */
            pos1 = AuxExonPtr->Position1;
            pos2 = AuxExonPtr->Position2;

            CheckBoundaries(&pos1, &pos2, length);

            if (pos2 - pos1 + 1 > MAXAA * LENGTHCODON) {
                printError
                    ("Not enough memory to hold exons : change MAXAA parameter");
            }

            lengthCDS += pos2 - pos1 + 1;

            if (lengthCDS > MAXAA * LENGTHCODON) {
                printError
                    ("Not enough memory to hold cds : change MAXAA parameter");
            }

            /* Printing all exons: DEFAULT */
            if (!(ALLE)) {
                printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
            }

            /* Joining all the CDS fragments */
            if (CDS) {
                if ((AuxExonPtr)->Strand == '-') {
                    ReverseSubSequence(pos1 - 1, pos2 - 1, Sequence, AuxExon);
                    AuxExon[pos2 - pos1 + 1] = '\0';
                    /* Joining the exons to build the complete CDS in reverse strand */
                    strcpy(AuxExon_tmp, AuxExon);
                    strcat(AuxExon_tmp, AuxGene);
                    strcpy(AuxGene, AuxExon_tmp);
                }
                else {
                    strncpy(AuxExon, Sequence + pos1 - 1, pos2 - pos1 + 1);
                    AuxExon[pos2 - pos1 + 1] = '\0';
                    /* Joining the exons to build the complete CDS in forward strand */
                    strcat(AuxGene, AuxExon);
                }
            }

            /* Printing sites or exons selected by type */

            /* Reverse strand */
            if (AuxExonPtr->Strand == '-') {
                /* First Exons */
                if (!strcmp(AuxExonPtr->Type, SFIRST)) {
                    /* Printing only first */
                    if (EFP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* START Site */
                    if (ALLS || SFP) {
                        pos1_site = pos2 - DOWN_START;
                        pos2_site = pos2 + UP_START;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTART);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }

                    /* DONOR Site */
                    if (ALLS || SDP) {
                        pos1_site = pos1 - DOWN_DONOR;
                        pos2_site = pos1 + UP_DONOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SDONOR);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }
                }

                /* Internal Exons */
                else if (!strcmp(AuxExonPtr->Type, SINTERNAL)) {
                    /* Printing only internal */
                    if (EIP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* ACCEPTOR Site */
                    if (ALLS || SAP) {
                        pos1_site = pos2 - DOWN_ACCEPTOR;
                        pos2_site = pos2 + UP_ACCEPTOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SACCEPTOR);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }

                    /* DONOR Site */
                    if (ALLS || SDP) {
                        pos1_site = pos1 - DOWN_DONOR;
                        pos2_site = pos1 + UP_DONOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SDONOR);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }
                }
                /* Terminal Exons */
                else if (!strcmp(AuxExonPtr->Type, STERMINAL)) {
                    /* Printing only terminal */
                    if (ETP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* Acceptor Site */
                    if (ALLS || SAP) {
                        pos1_site = pos2 - DOWN_ACCEPTOR;
                        pos2_site = pos2 + UP_ACCEPTOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SACCEPTOR);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }

                    /* Stop Site */
                    if (ALLS || STP) {
                        pos1_site = pos1 - DOWN_STOP;
                        pos2_site = pos1 + UP_STOP;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTOP);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }
                }

                /* Single Exons */
                else if (!strcmp(AuxExonPtr->Type, SSINGLE)) {
                    /* Printing only Single */
                    if (ESP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* START Site */
                    if (ALLS || SFP) {
                        pos1_site = pos2 - DOWN_START;
                        pos2_site = pos2 + UP_START;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTART);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }

                    /* STOP Site */
                    if (ALLS || STP) {
                        pos1_site = pos1 - DOWN_STOP;
                        pos2_site = pos1 + UP_STOP;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTOP);
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }
                }
                /* Not defined exon types */
                else {
                    /* First Site */
                    if (ALLS || SFP) {
                        pos1_site = pos1 - DOWN_SITE2;
                        pos2_site = pos1 + UP_SITE2;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, "Site_2");
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }

                    /* Second Site */
                    if (ALLS || SFP) {
                        pos1_site = pos2 - DOWN_SITE1;
                        pos2_site = pos2 + UP_SITE1;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, "Site_1");
                        printFasta(pos1_site, pos2_site, '-', Sequence);
                    }
                }
            }
            /* Forward Strand */
            else {
                /* First Exons */
                if (!strcmp(AuxExonPtr->Type, SFIRST)) {
                    /* Printing only first */
                    if (EFP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* START Site */
                    if (ALLS || SFP) {
                        pos1_site = pos1 - UP_START;
                        pos2_site = pos1 + DOWN_START;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTART);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }

                    /* DONOR Site */
                    if (ALLS || SDP) {
                        pos1_site = pos2 - UP_DONOR;
                        pos2_site = pos2 + DOWN_DONOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SDONOR);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }
                }

                /* Internal Exons */
                else if (!strcmp(AuxExonPtr->Type, SINTERNAL)) {
                    /* Printing only Internal exons */
                    if (EIP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* ACCEPTOR Site */
                    if (ALLS || SAP) {
                        pos1_site = pos1 - UP_ACCEPTOR;
                        pos2_site = pos1 + DOWN_ACCEPTOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SACCEPTOR);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }

                    /* DONOR Site */
                    if (ALLS || SDP) {
                        pos1_site = pos2 - UP_DONOR;
                        pos2_site = pos2 + DOWN_DONOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SDONOR);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }
                }
                /* Terminal Exons */
                else if (!strcmp(AuxExonPtr->Type, STERMINAL)) {
                    /* Printing only terminal */
                    if (ETP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* ACCEPTOR Site */
                    if (ALLS || SAP) {
                        pos1_site = pos1 - UP_ACCEPTOR;
                        pos2_site = pos1 + DOWN_ACCEPTOR;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SACCEPTOR);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }

                    /* STOP Site */
                    if (ALLS || STP) {
                        pos1_site = pos2 - UP_STOP;
                        pos2_site = pos2 + DOWN_STOP;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTOP);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }
                }

                /* Single Exons */
                else if (!strcmp(AuxExonPtr->Type, SSINGLE)) {
                    /* Printing only first */
                    if (ESP) {
                        printHeadExon(AuxExonPtr->Group, long_cc, Locus);
                        printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
                    }

                    /* START Site */
                    if (ALLS || SFP) {
                        pos1_site = pos1 - UP_START;
                        pos2_site = pos1 + DOWN_START;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTART);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }

                    /* STOP Site */
                    if (ALLS || STP) {
                        pos1_site = pos2 - UP_STOP;
                        pos2_site = pos2 + DOWN_STOP;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, SSTOP);
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }
                }
                else {
                    /* Print sites if Exon is not defined */
                    /* First Site */
                    if (ALLS) {
                        pos1_site = pos1 - UP_SITE1;
                        pos2_site = pos1 + DOWN_SITE1;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, "Site_1");
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }

                    /* Second Site */
                    if (ALLS) {
                        pos1_site = pos2 - UP_SITE2;
                        pos2_site = pos2 + DOWN_SITE2;
                        CheckBoundaries(&pos1_site, &pos2_site, length);
                        printHeadSite(AuxExonPtr->Group, long_cc, Locus, "Site_2");
                        printFasta(pos1_site, pos2_site, '+', Sequence);
                    }
                }
            }

        }

        /* Printing complete CDS */
        if (CDS) {
            printFastaGene(AuxExonPtr->Group, Locus, AuxGene);
            /* Initialize AuxGene array */
            AuxGene[0] = '\0';

        }

        /* Printing Transcripts */
        if (TRN) {
            /* Defining exon position1 and position2 */
            pos1 = ((genes + long_gg)->First->Strand == '+') ?
                   (genes + long_gg)->First->Position1 - UP_TRANSCRIP : (genes
                                                                         + long_gg)->First->
                   Position1 - DOWN_TRANSCRIP;
            pos2
                = ((genes + long_gg)->First->Strand ==
                   '+') ? (genes + long_gg)->Last->Position2 + DOWN_TRANSCRIP : (genes
                                                                                 + long_gg)->
                  Last->Position2 + UP_TRANSCRIP;

            CheckBoundaries(&pos1, &pos2, length);

            /* Extract the squence and stored them in AuxTrans */
            printHeadTranscript((genes + long_gg)->First->Group, Locus);
            printFasta(pos1, pos2, (genes + long_gg)->First->Strand, Sequence);
        }
    }

    /* 20150917 dk  mem leak fix1 */

    if (AuxExon_tmp != NULL) {
        free(AuxExon_tmp);
        AuxExon_tmp = NULL;
    }

    if (AuxExon != NULL) {
        free(AuxExon);
        AuxExon = NULL;
    }

    if (AuxGene != NULL) {
        free(AuxGene);
        AuxGene = NULL;
    }

    if (dd_dict != NULL) {
        free(dd_dict);
        dd_dict = NULL;
    }

    if (genes != NULL) {
        free(genes );
        genes = NULL;
    }

    if (exons != NULL) {
        free(exons );
        exons = NULL;
    }

    if (Sequence != NULL) {
        free(Sequence );
        Sequence = NULL;
    }

    /*
if (AuxExonPtr_tmp != NULL)
    {
      free (AuxExonPtr_tmp);
      AuxExonPtr_tmp = NULL;
    }

if (AuxExonPtr != NULL)
    {
      free (AuxExonPtr);
      AuxExonPtr = NULL;
    }

dk:  bad pointers -> unused above? here.... */

    return (0);
}
