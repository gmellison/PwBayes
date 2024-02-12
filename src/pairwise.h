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


typedef struct pairwisedists 
    {
    int           nTaxa;
    int           nPairs;
    MrBFlt        *dists;    /*  pointers to distances between taxa */
/* int        *doubletCounts;   */   
    int           *clusters; /*  not used  */
    int           *index;
    } PairwiseDists;

typedef int (*TiProbFxn_Pw)(PairwiseDists *, int, int);

/* 
extern int    numPairs;
extern int    numTrips;
extern int    *pairwiseCounts;
extern int    **tripleCounts;
extern int    defPairwise;
extern int    defTriples;
extern int    usePairwise;
extern int    useTriples;
 */

/* extern int    **tripletCounts;  */

void CountFreqs(void);
int  CountDoublets(int nPairs);
int  CompressDoubletData(void);
/* int  CountTriples(void); */
int  triplePos(int i, int j, int k);

int DoEstQPairwise(void);
int DoPairwiseLogLike(void);
int DoTripletLogLike(void);

int  DoPwSetParm(char *parmName, char *tkn); // param set -- see param list at end of command.c
void PrintPairwiseDists(PairwiseDists *pd);

PairwiseDists* AllocatePairwiseDists(void);
int FreePairwiseDists(PairwiseDists* pd);

void InitPairwiseDists(Tree *tree, PairwiseDists  *pd);
void InitPairwiseDistsPolyTree(PolyTree *tree, PairwiseDists *pd);

int TiProbsPairwise_JukesCantor (int division, int chain);
int DoubletProbs_JukesCantor(int division, int chain);
int TiProbsPairwise_Gen (int division, int chain);
int DoubletProbs_Gen(int division, int chain);
int CalcPairwiseDists_ReverseDownpass(Tree *t, int division, int chain);
int CountPairwise(void);
int Likelihood_Pairwise(int division, int chain, MrBFlt *lnL);

int CountTriplets(void);

int FreePairwise(void);
int FreeTriples(void);

int CalcTripletCnDists(int division, int chain);
int TiProbsTriplet_JukesCantor (int division, int chain);
int TripletProbs_JukesCantor (int division, int chain);
int Likelihood_Triples(int division, int chain, MrBFlt *lnL);

MrBFlt LogLikePairwise(int chain);
MrBFlt LogLikeTriplet(int chain);
MrBFlt LogLikeTriplet_Alpha(int chain);

#endif /* end of include guard PAIRWISE_H */

