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

#ifndef __PAIRWISE_H__
#define __PAIRWISE_H__

int TiProbsPairwise_JukesCantor (int division, int chain);
int DoubletProbs_JukesCantor(int division, int chain);
int TiProbsPairwise_Gen (int division, int chain);
int DoubletProbs_Gen(int division, int chain);
int CalcPairwiseDists_ReverseDownpass(Tree *t, int division, int chain);
int CountPairwise(int division);
int CalcPairwiseWeights(int chain);
int CalcPairwiseWeights_GTR(int chain);
int Likelihood_Pairwise(int division, int chain, MrBFlt *lnL);

int InitPairwise(void);
int FreePairwise(int numCurrentDivisions);

MrBFlt LogLikePairwise(int chain);

int PrepareHybridStep(int chain);
int PostHybridStep(int chain);

#endif /* end of include guard PAIRWISE_H */

