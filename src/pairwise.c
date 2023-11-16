/*
 * =====================================================================================
 *
 *       Filename:  pairwise.c
 *
 *    Description: G 
 *
 *        Version:  1.0
 *        Created:  11/16/2023 12:59:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
/**
 * @author      : greg (greg@$HOSTNAME)
 * @file        : pairwise
 * @created     : Thursday Nov 16, 2023 12:59:45 EST
 */

#include "pairwise.h"
#include "bayes.h"
#include "command.h"
#include "mbbeagle.h"
#include "model.h"
#include "mcmc.h"
#include "sumpt.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

// local prototypea:
int toIdx(int x);

void CountFreqs() {

    int s,i,j,nucIdx;
    int nucCounts[4] = {0};

    for (s=0;s<numChar;s++) {
        for (i=0;i<numTaxa;i++) {
            if (matrix[pos(i,s,numChar)]==GAP)
                continue;

            nucIdx=toIdx(matrix[pos(i,s,numChar)]);
            nucCounts[nucIdx]++;
        }
    }

    for (j=0;j<4;j++) {
        nucFreqs[j]=((MrBFlt)nucCounts[j])/(numChar*numTaxa);
    }

    defFreqs=YES;
}


void CountDoublets() {

    int s,i,j,nucIdx;
    int id1, id2;
    
    for (s=0;s<numChar;s++) {
        for (i=0;i<(numTaxa-1);i++) {
            if (matrix[pos(i,s,numChar)]==GAP)
                continue;

            id1=toIdx(matrix[pos(i,s,numChar)]);
            for (j=(i+1);j<numTaxa;j++) {
                 if (matrix[pos(j,s,numChar)]==GAP)
                    continue;
                    
                 id2=toIdx(matrix[pos(j,s,numChar)]);
                 doubletCounts[id1][id2]++;
            }
        }
    }
    for (i=0;i<4;i++) {
        for (j=1;j<4;j++) {
            doubletFreqs[i][j] = (doubletCounts[i][j] + doubletCounts[j][i]) / (2.0*numTaxa*numChar);
        }
    }
 
    defDoublets=YES;
}

int toIdx(int x){
    if (x == 1) return 0;
    else if (x == 2) return 1;
    else if (x == 4) return 2;
    else if (x == 8) return 3;
    else if (x == GAP) return 5;
    else {
        MrBayesPrint("Problem in to Idx (maybe bad input: x=%d?)\n", x);
        return ERROR;
    }
}

int triplePos(int i, int j, int k) {
    return i*16 +j*4 + k;
}


void CountTriples() {

    int s, seqi, seqj, seqk;
    int i,j,k;

    for (s=0;s<numChar;s++) {
        for (seqi=0;seqi<(numTaxa-2);seqi++) {
            if (matrix[pos(seqi,s,numChar)]==GAP)
                continue;
            i=toIdx(matrix[pos(seqi,s,numChar)]);

            for (seqj=(seqi+1);seqj<(numTaxa-1);seqj++) {
                 if (matrix[pos(seqj,s,numChar)]==GAP)
                    continue;
                 j = toIdx(matrix[pos(seqj,s,numChar)]);

                for (seqk=(seqj+1);seqk<numTaxa;seqk++) {
                     if (matrix[pos(seqk,s,numChar)]==GAP)
                        continue;
                     k = toIdx(matrix[pos(seqk,s,numChar)]);
                     tripleCounts[triplePos(i,j,k)]++;
                }
            }
        }
    }
    defTriples=YES;
}
