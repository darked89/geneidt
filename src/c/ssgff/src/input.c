/*************************************************************************
*                                                                        *
*   Module: input                                                        *
*                                                                        *
*   Formatted input of SSgff  (GFF and fasta)                            *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Genis PARRA FARRE                             *
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

extern int VRB;

FILE *OpenFile(char *FileName) {
    FILE *file;
    char mess[MAXSTRING];

    sprintf(mess, "File %s cannot be open for read", FileName);

    if ((file = fopen(FileName, "r")) == NULL) {
        printError(mess);
    }

    return(file);
}

long ReadSequence(FILE *seqfile, char *Locus, char *Sequence) {
    long pos;
    int  res;
    char cAux;

    /* line has the locus of the current sequence */
    res = fscanf(seqfile, ">%s", Locus);

    /* Jumping until \n of the first fasta line */
    res = fscanf(seqfile, "%c", &cAux);

    while (cAux != '\n') {
        res = fscanf(seqfile, "%c", &cAux);
    }

    /* fasta format = "atcgata...atta\n" */
    pos = 0;
    res = fscanf(seqfile, "%s\n", Sequence);

    /* Only one sequence can be read */
    while (res != EOF) {
        pos = pos + strlen(Sequence + pos);
        res = fscanf(seqfile, "%s\n", Sequence + pos);

        if ( !(pos % 10000000) ) {
            fprintf(stderr, "...%ld bp\n", pos);
        }

        if (pos >= SEQLENGTH) {
            printError("Not enough memory: change SEQLENGTH parameter ");
        }

    }

    /* End of sequence */
    pos = pos + strlen(Sequence + pos);
    fclose(seqfile);
    return(pos);
}

long ReadExonsGFF(FILE *file, exonGFF *exons, pack_gene *genes, dict *dd_dict, long *numgen, long Length) {
    char line[MAXLINE];
    char cc;
    long lastPos1;
    char mess[MAXSTRING];

    char *gff_column_1;
    char *gff_column_2;
    char *gff_column_3;
    char *gff_column_4;
    char *gff_column_5;
    char *gff_column_6;
    char *gff_column_7;
    char *gff_column_8;
    char *gff_column_9;

    long ii, gg;

    /* Coments: line begins with # */
    /* gff format = "Name  Source  Type  Begin  End  Score  Strand  Frame group */
    ii       = 0;
    gg       = 0;

    lastPos1 = -INFI;

    /* Skip comment lines */
    while (fgets(line, MAXLINE, file) != NULL) {

        if (line[0] == '#' || line[0] == '\n') {
            /* Skip this line */
            printMess("Skipping comment line");
        }
        else {
            /* For each line extract the features (GFF format) */
            gff_column_1 = (char *) strtok(line, "\t");
            gff_column_2 = (char *) strtok(NULL, "\t");
            gff_column_3 = (char *) strtok(NULL, "\t");
            gff_column_4 = (char *) strtok(NULL, "\t");
            gff_column_5 = (char *) strtok(NULL, "\t");
            gff_column_6 = (char *) strtok(NULL, "\t");
            gff_column_7 = (char *) strtok(NULL, "\t");
            gff_column_8 = (char *) strtok(NULL, "\t");
            gff_column_9 = (char *) strtok(NULL, "\n");

            if (gff_column_1 == NULL || gff_column_2 == NULL || gff_column_3 == NULL
                || gff_column_4 == NULL || gff_column_5 == NULL || gff_column_6 == NULL
                || gff_column_7 == NULL || gff_column_8 == NULL || gff_column_9 == NULL) {
                sprintf(mess, "Bad format: Exon GFF %ld\n", ii);
                printError(mess);
            }

            /* Exon Strand */
            if (sscanf(gff_column_7, "%c", &(exons[ii].Strand)) != 1) {
                sprintf(mess, "Bad format Strand: Exon %ld\n", ii);
                printError(mess);
            }

            /* 1/2. Sequence and Source not used */
            /* 3. Exon Type */
            if (sscanf(gff_column_3, "%s", exons[ii].Type) != 1) {
                sprintf(mess, "Bad format Type: Exon %ld\n", ii);
                printError(mess);
            }

            /* 4. Left position */
            if (sscanf(gff_column_4, "%ld", &(exons[ii].Position1)) != 1) {
                sprintf(mess, "Bad format Position1: Exon %ld\n", ii);
                printError(mess);
            }

            /* 4.b. Filename sort by position */
            if (exons[ii].Position1 < lastPos1) {
                sprintf(mess, "Bad position(not sorted): Exon %ld\n", ii);
                printError(mess);
            }

            lastPos1 = exons[ii].Position1;

            /* 5. Right position */
            if (sscanf(gff_column_5, "%ld", &(exons[ii].Position2)) != 1) {
                sprintf(mess, "Bad format Position2: Exon %ld\n", ii);
                printError(mess);
            }

            /* 5.b Boundaries position */
            if ((exons[ii].Position1) < 0 || (exons[ii].Position1) > Length
                || (exons[ii].Position2) < 0 || (exons[ii].Position2) > Length) {
                sprintf(mess, "Coordinates out of range: Exon %ld\n", ii);
                printError(mess);
            }

            if ((exons[ii].Position1) > (exons[ii].Position2)) {
                sprintf(mess, "Bad format positions : Exon %ld\n", ii);
                printError(mess);
            }

            /* 6. Score = '.' or float */
            if (sscanf(gff_column_6, "%lf", &(exons[ii].Score)) != 1) {
                if ((sscanf(gff_column_6, "%c", &cc) != 1) || (cc != '.')) {
                    sprintf(mess, "Bad format Score: Exon %ld\n", ii);
                    printError(mess);
                }

                exons[ii].Score = NOSCORE;
            }

            /* 7. Strand was done before */

            /* 8. Frame = '.' or integer */
            if (sscanf(gff_column_8, "%hd", &(exons[ii].Frame)) != 1) {
                if ((sscanf(gff_column_8, "%c", &cc) != 1) || (cc != '.')) {
                    sprintf(mess, "Bad format Frame: Exon %ld\n", ii);
                    printError(mess);
                }

                exons[ii].Frame = NOFRAME;
            }

            /* 9. Group */
            if (sscanf(gff_column_9, "%s", exons[ii].Group) != 1) {
                sprintf(mess, "Bad format Group: Exon %ld\n", ii);
                printError(mess);
            }

            /* Filling Gene structure */

            /* Assign a integer to the group name */
            gg = getkeyDict(dd_dict, (exons + ii)->Group);

            /* If it does not exists assign a new integer to the gene name */
            if (gg == NOTFOUND) {
                setkeyDict(dd_dict, (exons + ii)->Group);
                gg                     = getkeyDict(dd_dict, (exons + ii)->Group);
                /* Initializing a new gene */
                (genes + gg)->First    = (exons + ii);
                (genes + gg)->numExons = 0;
                (genes + gg)->Last     = (exons + ii);
                /* Increment the number of genes */
                (*numgen)++;
            }
            else {
                /* Assigning the to the previous exon of
                 this gene (genes+gg) the NextExon variable
                 pointing to the current exon (exons+ii)*/
                (genes + gg)->Last->NextExon = exons + ii;
                /* Increment numExons and assign the new Last exon*/
                (genes + gg)->numExons++;
                (genes + gg)->Last = (exons + ii);
            }

            /* printf ("%d %d %d %d \n",*numgen, (genes+gg)->numExons,
                 (genes+gg)->First->Position1,(genes+gg)->Last->Position1); */

            /* Process features from current exon */
            if (ii > NUMEXONS) {
                printError("Too many exons: Change NUMEXONS parameter");
            }

            ii++;

        }
    }
    fclose(file);
    return(ii);
}

