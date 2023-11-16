/*
 * =====================================================================================
 *
 *       Filename:  pairwise.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/16/2023 13:00:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
/**
 * @author      : greg (greg@$HOSTNAME)
 * @file        : pairwise
 * @created     : Thursday Nov 16, 2023 13:00:26 EST
 */

#ifndef PAIRWISE_H

#define PAIRWISE_H

#include "bayes.h"
#include "command.h"
#include "mbbeagle.h"
#include "model.h"
#include "mcmc.h"
#include "sumpt.h"
#include "utils.h"

MrBFlt nucFreqs[4];
int doubletCounts[4][4] = {{0}};
MrBFlt doubletFreqs[4][4] = {{0.0}};
int tripleCounts[64] = {0};

int defFreqs=NO;
int defDoublets=NO;
int defTriples=NO;

void CountFreqs();
void CountDoublets();
void CountTriples();
int  triplePos(int i, int j, int k);

#endif /* end of include guard PAIRWISE_H */

