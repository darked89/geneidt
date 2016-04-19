/*************************************************************************
*                                                                        *
*   Module: Dictionary                                                   *
*                                                                        *
*   Implementation of a look-up table by using hash tables               *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
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

/*  $Id: Dictionary.c,v 1.1 2000/07/05 08:12:02 eblanco Exp $  */

#include "ssgff.h"

/* Initializing the dictionary: hash table and counter of keys */
void resetDict(dict *my_hash) {
    int ii;

    for (ii = 0; ii < MAXENTRY; ii++) {
        my_hash->T[ii] = NULL;
    }
    my_hash->nextFree = 0;
}

/* Hash Function:: String -> Integer between 0..MAXENTRY-1 */
int hash_func(char s[]) {
    unsigned int ii;
    int          total;

    total = 0;

    for (ii = 0; ii < strlen(s); ii++) {
        total = (ii + 1) * s[ii] + total;
    }
    total = total % MAXENTRY;

    return(total);
}

/* Assign a number-key to the new word and store it */
int setkeyDict(dict *my_hash, char s[]) {
    int  ii;
    int  key;
    node *p_node;
    node *n_node;

    /* If this word exists at the dictionary don't insert */
    key = getkeyDict(my_hash, s);

    if (key == NOTFOUND) {
        ii = hash_func(s);

        /* Alloc the new word */
        if ((n_node = (node *) malloc(sizeof(node))) == NULL) {
            printError("Not enough memory: dictionary word");
        }

        /* Filling the node */
        strcpy(n_node->s, s);
        n_node->key = my_hash->nextFree++;

        if (my_hash->T[ii] == NULL) {
            n_node->next   = NULL;
            my_hash->T[ii] = n_node;
        }
        else {
            /* There are more nodes in this position: Colission */
            p_node         = my_hash->T[ii];
            /* Insert at the begining of the list */
            my_hash->T[ii] = n_node;
            n_node->next   = p_node;
        }

        key = n_node->key;
    }

    return(key);
}

/* Returns the key for the word request; NOTFOUND is Not found */
int getkeyDict(dict *my_hash, char s[]) {
    int  ii;
    int  key;
    int  found = 0;

    node *p_node;

    /* Computing hash function */
    ii = hash_func(s);

    /* Empty list means not found */
    if (my_hash->T[ii] == NULL) {
        key = NOTFOUND;
    }
    else {
        /* There are more nodes in this position: run the list */
        p_node = my_hash->T[ii];

        /* Searching until the first position not used */
        while ( p_node != NULL && !found ) {
            /* Same hash value: compare to see if it is the same string */
            if (!strcmp(s, p_node->s)) {
                found = 1;
                key   = p_node->key;
            }

            p_node = p_node->next;
        }

        if (!found) {
            key = NOTFOUND;
        }
    }

    return(key);
}

/* Shows the dictionary */
void showDict(dict *my_hash) {
    int  ii;
    node *p_node;

    printf("Dictionary: \n\n");
    for (ii = 0; ii < MAXENTRY; ii++) {
        if (my_hash->T[ii] != NULL) {
            /* There are more nodes in this position */
            p_node = my_hash->T[ii];

            /* Searching the first position free */
            while ( p_node != NULL ) {
                printf("%-20s | \t\t\t %d\n", p_node->s, p_node->key);
                p_node = p_node->next;

            }
        }
    }
}

/* Free memory of hash nodes (sinonimous) */
void freeNodes(pnode node) {
    if (node == NULL) {
    }
    else {
        freeNodes(node->next);
        free(node);
    }
}

/* Free memory of the whole dictionary */
void freeDict(dict *my_hash) {
    int ii;

    /* free all of the words in the dictionary */
    for (ii = 0; ii < MAXENTRY; ii++) {
        freeNodes(my_hash->T[ii]);
        my_hash->T[ii] = NULL;
    }

    /* free the dictionary */
    free(my_hash);
    my_hash = NULL;
}

/* Binding the amino acid (key) to the new codon (word) */

void setAADict(dict *my_hash, char sCodon[], char aA) {
    int  ii;
    node *p_node;
    node *n_node;

    ii = hash_func(sCodon);

    /* Allocating the new word */
    if ((n_node = (node *) malloc(sizeof(node))) == NULL) {
        printError("Not enough memory: AA dictionary node");
    }

    /* printf("n_node is at memory location %p:\n", n_node); */

    /* Filling the node */
    /* printf("before strcpy  %s:\n", sCodon); */

    strcpy(n_node->s, sCodon);
    /* dk 20150917a  */
    /*n_node->key = my_hash->nextFree++; */
    /* end dk 20150917a  */
    n_node->key = aA;

    if (my_hash->T[ii] == NULL) {
        n_node->next   = NULL;
        my_hash->T[ii] = n_node;
    }
    else {
        /* There are more nodes in this position: Colission */
        p_node         = my_hash->T[ii];
        /* Insert at the begining of the list */
        my_hash->T[ii] = n_node;
        n_node->next   = p_node;
    }
}

/* Returns the amino acid for the input codon; '?' is Not found */
char getAADict(dict *my_hash, char s[]) {
    int  ii;
    int  aa;
    int  found = 0;

    node *p_node;

    aa = '?'; /* dk testing */

    ii = hash_func(s);

    if (my_hash->T[ii] == NULL) {
        aa = '?';
    }
    else {
        /* There are more nodes in this position */
        p_node = my_hash->T[ii];

        /* Searching until the first position free */
        while ( p_node != NULL && !found ) {
            if (!strcmp(s, p_node->s)) {
                found = 1;
                aa    = p_node->key;
            }

            p_node = p_node->next;
        }

        if (!found) {
            aa = '?';
        }
    }

    return(aa);
}
