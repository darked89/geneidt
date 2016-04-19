/*************************************************************************
*                                                                        *
*   Module: output                                                       *
*                                                                        *
*   Formatted output of SSgff                            )               *
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

/* Printing error messages */
void printError(char *s) {
    fprintf(stderr, "Error: %s\n", s);
    exit(1);
}

/* Printing messages (information) */
void printMess(char *s) {
    if (VRB) {
        fprintf(stderr, "> %s\n", s);
    }
}

void printFasta(long pos1, long pos2, char strand, char *seq) {
    long jj;
    long ii = 1;

    /* For reverse strand */
    if (strand == '-') {
        for (jj = pos2 - 1; jj >= pos1 - 1; jj--, ii++) {
            printf("%c", complement(seq[jj]));

            if (!(ii % 60)) {
                printf("\n");
            }
        }

    }
    else {
        /* Forward  strand */
        for (jj = pos1 - 1; jj <= pos2 - 1; jj++, ii++) {
            printf("%c", seq[jj]);

            if (!(ii % 60)) {
                printf("\n");
            }
        }

    }

    if ((ii - 1) % 60) {
        printf("\n");
    }

}

void printHeadSite(char *Group, int numexon, char *nameseq, char *Type) {

    printf(">%s.%d:%s %s\n", Group, numexon + 1, nameseq, Type);

}

void printHeadExon(char *Group, int numexon, char *nameseq) {

    printf(">%s.%d:%s exon %d \n", Group, numexon + 1, nameseq, numexon + 1);

}

void printFastaGene(char *Group, char *nameseq, char *AuxGene) {
    unsigned int jj;

    printf(">%s:%s CDS\n", Group, nameseq);

    for (jj = 0; jj < strlen(AuxGene); jj++) {
        printf("%c", AuxGene[jj]);

        if (!((jj + 1) % 60)) {
            printf("\n");
        }
    }

    if (jj % 60) {
        printf("\n");
    }
}

void printHeadTranscript(char *Group, char *nameseq) {

    printf(">%s:%s Primary Transcript \n", Group, nameseq);

}

void printHeadIntron(char *Group, int numexon, char *nameseq) {

    printf(">%s.i%d:%s Intron between exon %d and %d \n", Group, numexon, nameseq, numexon, numexon + 1);

}

