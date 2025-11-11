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
#include <math.h>
/**
 * @author      : greg (greg@$HOSTNAME)
 * @file        : pairwise
 * @created     : Thursday Nov 16, 2023 12:59:45 EST
 */

#include "bayes.h"
#include "pairwise.h"
#include "utils.h"
#include "command.h"
#include "mcmc.h"
#include "model.h"
#include "likelihood.h"

#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#define LIKEPW_EPSILON        1.0e-300
#define LIKEPW_EPSILON        1.0e-300

extern MrBFlt *pwWeight ;

// global variables for pairwise likelihoods:
MrBFlt nucFreqs[4];
MrBFlt doubletFreqs[4][4];

MrBFlt relRates[6];
int    tempNumStatesPw;

int defFreqs = NO;
int defDoublets = NO;

int *pairwiseCounts;
int numPairs;

/* globals declared here */
int defPairwise=NO;
int allocPairwise=NO;


// local prototypes:
int              toIdx(int x);
int              pairIdx(int i, int j, int n);

/*  *****************
 *
 *
 *  Utility Functions
 *
 *
 *  *****************  */

#define  dIdx( i,j, dim_j ) i*dim_j + j
#define  tIdx( k, i, j, dim_i, dim_j ) k*dim_i*dim_j + i*dim_j + j

int pairIdx(int i, int j, int n) {
    int k;

    if (i == j || i < 0 || j < 0 || i >= n || j >= n) {
        MrBayesPrint("Error in pair indexing; i=%d, j=%d,n=%d \n", i,j,n);
        return(-1);
    }

    if (i > j) {
        k = i;
        i = j;
        j = k;
    } 

    return(n * (n-1)/2 - (n-i)*(n-i-1)/2 + (j-i-1));
}

/*  *****************
 *
 *
 *  Pw likelihood for Mrb MCMC
 *   
 *
 *  *****************  */

int FreePairwise(int numLocalChains) {
    int         i, j;
    ModelInfo   *m;

    /* free model variables for Gibbs gamma */
    for (i=0; i<numCurrentDivisions; i++)
        {
        m=&modelSettings[i];

        /* free pairwise dists and transition probs  */        
        if (m->usePairwise == NO)  
            continue;

        if (m->pwDists)
            {
            for (j=0; j<numLocalChains; j++)
                {
                if (m->pwDists[j])
                    free(m->pwDists[j]);
                }
            free(m->pwDists);
            }

        if (m->tiProbsPw)
            {
            for (j=0; j<m->numTiProbsPw; j++)
                if (m->tiProbsPw[j]) /*  chain pw probs */
                    free(m->tiProbsPw[j]);
            free(m->tiProbsPw);
            }

        if (m->doubletProbs)
            {
            for (j=0; j<m->numDoubletProbs; j++)
                if (m->doubletProbs[j]) /*  chain pw probs */
                    free(m->doubletProbs[j]);
            free(m->doubletProbs);
            }

        if (m->pwIndex)
            {
            for (j=0; j<numLocalChains; j++)
                {
                if (m->pwIndex[j]) /*  chain pw probs */
                    free(m->pwIndex[j]);
                }
            free(m->pwIndex);
            }

        }
    return (NO_ERROR);
}

int CalcPairwiseDists_ReverseDownpass(Tree *t, int division, int chain)
{
    int         i,j,a,d,k;
    TreeNode    *p;
    double      x;
    MrBFlt      *dists, **distsTemp;
    //int         *distShare;
    ModelInfo   *m;

    m = &modelSettings[division];

    dists = m->pwDists[chain];
    //distShare =m->pwDistsShare[chain];

    int numExtNodes = numLocalTaxa;

    /*  We'll calculate (in distsTemp) all node dists including internal nodes  */
    distsTemp=(MrBFlt**)malloc(t->nNodes*sizeof(MrBFlt*)); 
    for (k=0;k<t->nNodes;k++)
            distsTemp[k]=(MrBFlt*)malloc(t->nNodes*sizeof(MrBFlt));

    /*  make sure dists are init to 0  */
    for (i=0; i<t->nNodes; i++) 
        for (j=0; j<t->nNodes; j++)
                distsTemp[i][j]=0.0;

    /* loop over all nodes  */
    for (i=(t->nNodes)-2; i>=0; i--)
        {

        p = t->allDownPass[i];
        a = p->anc->index;   
        d = p->index;

        //MrBayesPrint("  Outer Loop: Node %d \n", i);

        /* find length */
        if (m->cppEvents != NULL)
            {
            x = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
            }
        else if (m->tk02BranchRates != NULL)
            {
            x = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
            }
        else if (m->wnBranchRates != NULL)
            {
            x = GetParamSubVals (m->wnBranchRates, chain, state[chain])[p->index];
            }
        else if (m->ilnBranchRates != NULL)
            {
            x = GetParamSubVals (m->ilnBranchRates, chain, state[chain])[p->index];
            }
        else if (m->igrBranchRates != NULL)
            {
            x = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
            }
        else if (m->mixedBrchRates != NULL)
            {
            x = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
            }
        else
            x = p->length;

        /*  start with distance from node to ancestor  */
        distsTemp[d][a] = x;
        distsTemp[a][d] = x;

        /* now revisit previously visited nodes, updating distances  */ 
        for (j=i+1; j<t->nNodes; j++)
            {
                //MrBayesPrint("  Inner Loop: Node %d \n", j);
                k=t->allDownPass[j]->index;
                if (k==a) 
                    continue;

                /*  dist from current node to prior node =  
                 *      dist from anc to prior node + dist from anc to current node  */
                distsTemp[d][k] = distsTemp[a][k] + x;
                distsTemp[k][d] = distsTemp[k][a] + x;

            }
        }

    /*  set the distances in the output array */
    for (i=0; i<(numExtNodes-1); i++) {
        for (j=i+1; j<numExtNodes; j++) {
            dists[pairIdx(i,j,numExtNodes)]=distsTemp[i][j];
        }
    }

    ///* check which pairwise dists are same   */
    //m->numUniqueDists[chain]=0;
    //for (i=0;i<(numLocalTaxa*(numLocalTaxa-1)/2);i++)
    //    {
    //    distShare[i]=i;
    //    for (j=0;j<i;j++) 
    //        {
    //        if (AreDoublesEqual(dists[i], dists[j], ETA))
    //            {
    //            distShare[i]=j;
    //            break;
    //            }
    //        }
    //    }
    //for (i=0;i<(numLocalTaxa*(numLocalTaxa-1)/2);i++)
    //    {
    //    if (distShare[i] == i) 
    //        m->numUniqueDists[chain]+=1;
    //    }
    ////MrBayesPrint("   uniquedists: %d \n", m->numUniqueDists[chain]);

    /*  free temp array and return pointer to taxa pairwise distances */
    for (k=0;k<t->nNodes;k++)
        free(distsTemp[k]);
    free(distsTemp);

    return(NO_ERROR);
}

// int CalcPairwiseDistsPolyTree(PolyTree *t, PairwiseDists *pd)
// {
//     int         i,j,a,d,k;
//     PolyNode    *p;
//     double      x;
// 
//     int numExtNodes = t->nNodes - t->nIntNodes;
// 
//     /*  We'll calculate (in distsTemp) all node dists including internal nodes  */
//     MrBFlt **distsTemp;
//     distsTemp=(MrBFlt**)malloc(t->nNodes*sizeof(MrBFlt*)); 
//     for (k=0;k<t->nNodes;k++)
//             distsTemp[k]=(MrBFlt*)malloc(t->nNodes*sizeof(MrBFlt));
// 
//     /*  make sure dists are init to 0  */
//     for (i=0; i<t->nNodes; i++) 
//         for (j=0; j<t->nNodes; j++)
//                 distsTemp[i][j]=0.0;
// 
//     /* loop over all nodes  */
//     for (i=1; i<t->nNodes; i++)
//         {
//         p = &t->nodes[i];
//         a = p->anc->index;   
//         d = p->index;
//         x = p->length;
// 
//         /*  start with distance from node to ancestor  */
//         distsTemp[d][a] = x;
//         distsTemp[a][d] = x;
// 
//         /* now revisit previously visited nodes, updating distances  */ 
//         for (j=i-1; j>=0; j--)
//             {
//                 k=(&t->nodes[j])->index;
//                 if (k==a) 
//                 { 
//                     continue;
//                 }
// 
//                 /*  dist from current node to prior node =  
//                  *      dist from anc to prior node + dist from anc to current node  */
//                 distsTemp[d][k] = distsTemp[a][k] + x;
//                 distsTemp[k][d] = distsTemp[k][a] + x;
// 
//             }
//         }
// 
//     /*  set the distances in the output array */
//     for (i=0; i<(numExtNodes-1); i++) {
//         for (j=i+1; j<numExtNodes; j++) {
//             pd->dists[pairIdx(i,j,pd->nTaxa)]=distsTemp[i][j];
//         }
//     }
// 
//     /*  free temp array and return pointer to taxa pairwise distances */
//     for (k=0;k<t->nNodes;k++)
//         free(distsTemp[k]);
//     free(distsTemp);
// 
//     return(NO_ERROR);
// }

int toIdx(int x){
    if (x == 1) return 0;
    else if (x == 2) return 1;
    else if (x == 4) return 2;
    else if (x == 8) return 3;
    else if (x == GAP) return -1;
    else if (x == MISSING) return -1;
    else {
        MrBayesPrint("Problem in to Idx (maybe bad input: x=%d?)\n", x);
        return ERROR;
    }
}

/*-----------------------------------------------------------------
|
|   TiProbs_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsPairwise_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    
    int         i, j, k, p, index;
    MrBFlt      *dists; 
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange;
    CLFlt       *tiP;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    /* find pw dists and transition probabilities */
    dists = m->pwDists[chain];

    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
   
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    for (p=0; p<m->numPairs; p++)
        {

        /* only calc tiprobs once -- points to prev calculated probs  */
        //if (m->pwDistsShare[chain][p] < p) 
        //    //continue;
        //    {
        //    m->tiProbsPw[m->pwIndex[chain][p]] = m->tiProbsPw[m->pwIndex[chain][m->pwDistsShare[chain][p]]];
        //    continue;
        //    }

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        length = dists[p];

        /* fill in values */
        index=0;
        for (k=0; k<m->numRateCats; k++)
            {

            t = length * baseRate * catRate[k];

            /* numerical errors will ensue if we allow very large or very small branch lengths,
            which might occur in relaxed clock models */
            if (t < TIME_MIN)
                {
                /* Fill in identity matrix */
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = 1.0;
                        else
                            tiP[index++] = 0.0;
                        }
                    }
                }
            else if (t > TIME_MAX)
                {
                /* Fill in stationary matrix */
                for (i=0; i<4; i++)
                    for (j=0; j<4; j++)
                        tiP[index++] = 0.25;
                }
            else
                {
                /* calculate probabilities */
                pChange   = (CLFlt) (0.25 - 0.25 * exp(-(4.0/3.0)*t));
                pNoChange = (CLFlt) (0.25 + 0.75 * exp(-(4.0/3.0)*t));

                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = pNoChange;
                        else
                            tiP[index++] = pChange;

                        }
                    }
                }
            }
        } /*  end loop over pairs */

        /*  Pw transition probabilities are done, now let's fill in the doublet probs: */
    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   TiProbsPairwise_Gen: update transition probabilities for 4by4
|       nucleotide model with nst == 6 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsPairwise_Gen (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    register int    i, j, k, n, p, s, index;
    MrBFlt          t, *catRate, baseRate, *eigenValues, *cijk, *bs,
                    EigValexp[64], sum, *ptr, theRate,
                    length;
    CLFlt           *tiP;
    ModelInfo       *m;
    MrBFlt          *dists; 

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    n = m->numModelStates; 

    /* find pw dists and transition probabilities */
    dists = m->pwDists[chain];

    /* get base rate */
    baseRate = GetRate (division, chain);
   
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
   
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    /* get eigenvalues and cijk pointers */
    eigenValues = m->cijks[m->cijkIndex[chain]];
    cijk        = eigenValues + (2 * n);

    for (p=0; p<m->numPairs; p++)
        {

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        length = dists[p];

        /* fill in values */
        index=0;
        for (k=0; k<m->numRateCats; k++)
            {

            t = length * baseRate * catRate[k];

            /* numerical errors will ensue if we allow very large or very small branch lengths,
            which might occur in relaxed clock models */
            if (t < TIME_MIN)
                {
                /* Fill in identity matrix */
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = 1.0;
                        else
                            tiP[index++] = 0.0;
                        }
                    }
                }
            else if (t > TIME_MAX)
                {
                /* Get base freq */
                bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

                /* Fill in stationary matrix */
                for (i=0; i<4; i++)
                    for (j=0; j<4; j++)
                        tiP[index++] = (CLFlt) bs[j];
                }
            else
                {
                /* We actually need to do some work... */
                for (s=0; s<4; s++)
                    EigValexp[s] =  exp(eigenValues[s] * t);

                ptr = cijk;
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        sum = 0.0;
                        for (s=0; s<4; s++)
                            sum += (*ptr++) * EigValexp[s];
                        tiP[index++] = (CLFlt) ((sum < 0.0) ? 0.0 : sum);
                        }
                    }
                }
            }
        } /*  end loop over pairs */
        /*  Pw transition probabilities are done, now let's fill in the doublet probs: */
    return (NO_ERROR);
}

/*-----------------------------------------------------------------
|
|   DoubletProbabilities_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/

int DoubletProbs_JukesCantor(int division, int chain)
{
    int         i, j, k, p, index, dpIdx;
    CLFlt       *tiP, *doubP;   
    ModelInfo   *m;
    MrBFlt      *bs;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    for (p=0; p<m->numPairs; p++)
        {

        /* only calc tiprobs once -- points to prev calculated probs  */
        //if (m->pwDistsShare[chain][p] < p) 
            //  continue;

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        /*  reset doublet probs to 0: */
        dpIdx = 0;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                doubP[dpIdx++] = 0.0;
        
        index=0; 
        for (k=0; k<m->numRateCats; k++) 
            {
            dpIdx=0;
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    doubP[dpIdx++] += 0.25 * tiP[index++]/((MrBFlt)m->numRateCats);
            }
        }
    return(NO_ERROR);
}

int DoubletProbs_Gen(int division, int chain)
{

    int         i, j, k, p, index, dpIdx;
    CLFlt       *tiP, *doubP;   
    ModelInfo   *m;
    MrBFlt      *bs;

    m = &modelSettings[division];
    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    for (p=0; p<m->numPairs; p++)
        {

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        /*  reset doublet probs to 0: */
        dpIdx = 0;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                doubP[dpIdx++] = 0.0;
        
        index=0; 
        for (k=0; k<m->numRateCats; k++) 
            {
            dpIdx=0;
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    doubP[dpIdx++] += bs[i] * tiP[index++]/((MrBFlt)m->numRateCats);
            }
        }
    return(NO_ERROR);
}



/*  
 *  Original doublet count setup. Used global pairwiseCounts. 
 *  Need to update to allow for MrB partitioning
 *  TODO: convert this into "InitPairwiseCounts"
 *  along with "InitPairwiseProbs" to run in the main MCMC method  
 *  */
int InitPairwise(void) {

    int             i,j,k,c,d,id1,id2;
    ModelInfo       *m;
    
    /*  For now, only inplemented for a single partition 
     *  of DNA type data */
    /*  
    MrBayesPrint("Counting Pairwise\n"); 
    if (mp->dataType!=DNA && mp->dataType != RNA)
        { 
        MrBayesPrint("%s Can't count pairwise site-patterns for non DNA data  \n", spacer);
        return(ERROR);
        }
    */
   if (memAllocs[ALLOC_PAIRWISE] == YES) 
       {
       free(pairwiseCounts);
       for (d=0; d<numCurrentDivisions; d++) 
           modelSettings[d].pwCounts=NULL;
       memAllocs[ALLOC_PAIRWISE] = NO;
       }


    for (d=0; d<numCurrentDivisions; d++) 
        {
        m=&modelSettings[d];
        if (!m->usePairwise)
            continue;

        // MrBayesPrint("Counting Pairwise\n");
        if (defMatrix == NO) 
            {
            MrBayesPrint("%s Matrix needs to be defined before counting doublets  \n", spacer);
            return(ERROR);
            }  
        
        //if (numDefinedPartitions > 1)    
        //    { 
        //    MrBayesPrint("%s Pairwise count likelihood only implemented for a single partition  \n", spacer);
        //    return(ERROR);
        //    }

        m->numPairs=numLocalTaxa*(numLocalTaxa-1)/2;
        /* first allocate pairwise doublet counts:  */
        m->pwCounts=(int*)SafeMalloc(m->numPairs * m->numStates * m->numStates * sizeof(int));
        if (!m->pwCounts)
            {
            MrBayesPrint("%s Problem allocating pairwise counts! \n", spacer);
            free(m->pwCounts);
            return(ERROR);
            }

        /* now count the doublets across taxa pairs */
        for (i=0; i<(numLocalTaxa-1); i++)
            {
            for (j=i+1; j<numLocalTaxa; j++)
                {
                k=pairIdx(i,j,numLocalTaxa);
                for (c=0;c<numChar;c++)
                    {
                    if (matrix[pos(i,c,numChar)]==GAP | matrix[pos(j,c,numChar)]==GAP | matrix[pos(i,c,numChar)]==MISSING | matrix[pos(j,c,numChar)]==MISSING )
                        continue;

                    if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != d+1)
                        continue;

                     /* nucleotides at position x of sequences i & j   */        
                    id1=toIdx(matrix[pos(i,c,numChar)]);
                    id2=toIdx(matrix[pos(j,c,numChar)]);

                    /* increment the count of that nucleotide pair, at pair k  */
                    m->pwCounts[tIdx(k,id1,id2,m->numStates,m->numStates)]++;

                    }
                }
            }
        }

    memAllocs[ALLOC_PAIRWISE] = YES;

    return (NO_ERROR);
}



/*
 * Utility function for resetting CI calculation flags. 
 */
int PrepareHybridStep(int chain)
{
    ModelInfo   *m;
    Tree        *tree;
    TreeNode    *p;
    int         i,d;

    for (d=0; d<numCurrentDivisions; d++)
    {
        m = &modelSettings[d];
        tree = GetTree(m->brlens, chain, state[chain]);
        m->usePairwise=NO;
        m->upDateCijk=YES;

        for (i=0; i<tree->nIntNodes; i++) 
            {
            p = tree->intDownPass[i];
            p->left->upDateTi=YES;
            p->right->upDateTi=YES;
            p->upDateCl=YES; 
            }
    }

    return (NO_ERROR);
}

int PostHybridStep(int chain)
{
    ModelInfo   *m;
    int         d;

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        m->usePairwise=YES;
        }

    return (NO_ERROR);
}
/*-----------------------------------------------------------------
|
|   Likelihood_Pairwise: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int Likelihood_Pairwise (int division, int chain, MrBFlt *lnL)
{
    /* calculate Jukes Cantor likelihood, using pw counts and tps. */
    
    int         i, j, idx, p, nijk;
    MrBFlt      like;
    CLFlt       *doubP, pijk;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    (*lnL)=0.0;

    for (p=0; p<m->numPairs; p++)
        {

        /* find transition probabilities */
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        idx=0;
        for (i=0; i<4; i++) 
            {
            for (j=0; j<4; j++) 
                {
                like = 0.0;
                nijk=m->pwCounts[tIdx(p,i,j,4,4)];
                pijk=doubP[idx++];

                if (nijk == 0)
                    like+=0;
                else 
                    like+=nijk*log(pijk);         
 
                /* check against LIKE_EPSILON (values close to zero are problematic) */
                if (pijk < LIKEPW_EPSILON)
                    {
                    (*lnL) = MRBFLT_NEG_MAX;
                    abortMove = YES;
                    return ERROR;
                    }

                if (m->usePwWeights > 0)
                    like=like*pwWeight[division];

                (*lnL)+=like;

                }
            }
        }

    return (NO_ERROR);
}

MrBFlt LogLikePairwise(int chain) 
{
    ModelInfo  *m;
    Tree       *tree;
    MrBFlt     chainLnLike;
    int d;

    chainLnLike = 0.0;
   
    for (d=0; d<numCurrentDivisions; d++) 
        {
        m = &modelSettings[d];
        tree = GetTree(m->brlens, chain, state[chain]);

        if (m->upDateCijk == YES)
            {
            if (UpDateCijk(d, chain)== ERROR)
                {
                (m->lnLike[2*chain+state[chain]]) = MRBFLT_NEG_MAX; /* effectively abort the move */
                return (MRBFLT_NEG_MAX);
                }
            m->upDateAll = YES;
            }

        CalcPairwiseDists_ReverseDownpass(tree,d,chain);
        m->PwTiProbs(d,chain);
        m->DoubletProbs(d,chain);

        m->PwLikelihood(d,chain,&(m->lnLike[2*chain+state[chain]]));
        chainLnLike+=m->lnLike[2*chain+state[chain]];
        }

    return(chainLnLike);
}

MrBFlt EstPwDist_JC(int p) 
{
    MrBFlt tau;
    tau = 0.0;   
    return(tau);
}


/* 
 * calculate the weights for adjusting composite likelihoods within the MH proposal. 
 * The weights depend on 1. the data and 2.  the max comp. like. estimate 
 *   also need the rate shape parameter to be given. 
 * */



//int PwWeights_JC(int d) 
//{
//    MrBFlt *pww;
//    ModelInfo *m;
//    
//    m = &modelSettings[d];
//    pww = m->pwLikeWeights;
//
//    /* first calculate the n_ii for the division */
//
//
//    return(0);
//}


/* ****************************
 * 
 *
 *  Calculations for weighting 
 *  Pw likelihood. 
 *
 *
   **************************** */

int CalcPairwiseWeights (int chain) {

    /*  setup the arays for pairwise counts,
     *  and populate the counts.
     * 
     *  this function will supercede the 'PairwiseCounts'
     *  function, since that just uses a global variable and 
     *  we want to be able to apply pw likelihood within partitions.
     *
     *  We'll also make and populate the counts for the data splits
     *  for estimating jacobians for weighting the pw likelihood. 
     *  */

    ModelInfo* m;

    int i,j,k,l,c,c1,c2,d,z;
    int cI,dI,pI,splitI;
    int nI, nIOverall;
    int *counts, count;
    int **countIndex;
    int overallPwIdx;
    int index, indexStep;
    int nSplits, nPairs, nStates;
    int *n10, *n11;
    MrBFlt *pwDists;
    MrBFlt baseRate;
    MrBFlt *catRate;
    MrBFlt theRate;
    MrBFlt rc;
    MrBFlt *p_10, *p_11, *p1_10, *p1_11; //, *p2_10, *p2_11; 
    int    nidx;
    MrBFlt **J, **H, **Hinv;
    MrBFlt **HiJ, *eigvals, *eigvalsc;
    int **niiIndex ;
    MrBFlt t1,t2;
    MrBFlt eterm;
    int numBranches;
    Tree *tree;
    int freeBitsets;
    int *tempPartitionPair;
    int countLen;

    int nLongsNeeded;
    int **PairBranch;
    int l1;
    int l2;
    int nRates;
    MrBFlt al;

    /*  first set up worker matrices for eigen computation */
    // set up matrices 
    MrBFlt **V, **Vinv;
    MrBComplex **Vc, **Vcinv; 

    MrBFlt* dw;
    int*    iw;
    int isComplex;
    MrBFlt eigsum;
    MrBFlt eigsum2;
    MrBFlt em;
    MrBFlt v;

    // loop over the partitions:
    for (d=0; d<numCurrentDivisions; d++)
        {

        m = &modelSettings[d];
        tree = GetTree(m->brlens, chain, state[chain]);
        if (m->usePwWeights==3)
            {
            pwWeight[d] = 2.0 / (numLocalTaxa * (numLocalTaxa - 1));
            MrBayesPrint("%s pwWeight: %f \n", spacer, pwWeight[d]);
            continue;
            }
                
        nStates = m->numModelStates;
        nSplits = m->numDataSplits;
        nPairs = m->numPairs;

        if (!tree->isRooted)
            numBranches=numLocalTaxa*2 - 3;
        else 
            numBranches=numLocalTaxa*2 - 2;

        overallPwIdx = m->numDataSplits;
        MrBayesPrint("%s Calculating PW weight for JC submod using method %d. \n", spacer, m->usePwWeights);

        V     = AllocateSquareDoubleMatrix(numBranches);
        Vinv  = AllocateSquareDoubleMatrix(numBranches);
             
        Vc    = AllocateSquareComplexMatrix(numBranches); 
        Vcinv = AllocateSquareComplexMatrix(numBranches);

        /*  * 
         *  Initialize necessary arrays:
         *  */
        /* initialize counts & index array */

        countLen = (nSplits+1) * nPairs * nStates * nStates;

        counts = (int*) SafeMalloc( (nSplits+1) * nPairs * nStates * nStates * sizeof(int));
        if (counts == NULL)
            return(ERROR);

        countIndex=(int**) SafeMalloc( (nSplits+1) * sizeof(int*));
        if (countIndex == NULL)
            return(ERROR);        

        for (i=0; i<(nSplits+1); i++)
            {
            countIndex[i]=(int*) SafeMalloc( nPairs * sizeof(int));
            if (countIndex[i] == NULL)
                return(ERROR);        
            }


        index=0;
        indexStep=nStates*nStates;
        for (i=0; i<(nSplits+1); i++)
            {
            for (j=0; j<nPairs; j++)
                {
                countIndex[i][j]=index;
                index+=indexStep;
                }
            }

        /*  init array for pw distances */
        pwDists = (MrBFlt*) SafeMalloc( nPairs * sizeof(MrBFlt));
        if (pwDists == NULL) 
            return (ERROR);

        n10 = (int*) SafeMalloc( (nSplits+1) * nPairs * sizeof(int));
        if (n10 == NULL) 
            return (ERROR);

        n11 = (int*) SafeMalloc( (nSplits+1) * nPairs * sizeof(int));
        if (n11 == NULL) 
            return (ERROR);

        niiIndex=(int**) SafeMalloc( (nSplits+1) * sizeof(int*));
        for (i=0; i<(nSplits+1); i++)
            niiIndex[i]=(int*) SafeMalloc( nPairs * sizeof(int));

        index=0;
        indexStep=1;
        for (splitI=0; splitI<(nSplits+1); splitI++)
            {
            for (j=0; j<nPairs; j++)
                {
                niiIndex[splitI][j]=index;
                index+=indexStep;
                }
            }

        /* init arrays for derivatives  */
        p_10 = (MrBFlt*) SafeMalloc( nPairs * sizeof(MrBFlt));
        if (p_10 == NULL)
            return(ERROR);
        p_11 = (MrBFlt*) SafeMalloc( nPairs * sizeof(MrBFlt));
        if (p_11 == NULL)
            return(ERROR);
        p1_10 = (MrBFlt*) SafeMalloc(nPairs * sizeof(MrBFlt));
        if (p1_10 == NULL)
            return(ERROR);
        p1_11 = (MrBFlt*) SafeMalloc(nPairs * sizeof(MrBFlt));
        if (p1_11 == NULL)
            return(ERROR);

        for (k=0;k<nPairs;k++)
            {
            p_10[k]  = 0.0;
            p_11[k]  = 0.0;
            p1_10[k] = 0.0;
            p1_11[k] = 0.0;
            }

        /*  init arrays for first and second derivs */
        MrBFlt  **D1L, **D1LP;
        D1L =  (MrBFlt**) SafeMalloc( nSplits * sizeof(MrBFlt*));
        if (!D1L)
            return (ERROR);

        for (i=0; i<nSplits; i++) 
            { 
            D1L[i]=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
            if (!D1L[i])
                return (ERROR);
            }

        D1LP =  (MrBFlt**) SafeMalloc( nSplits * sizeof(MrBFlt*));
        if (!D1LP)
            return (ERROR);

        for (i=0; i<nSplits; i++) 
            { 
            D1LP[i]=(MrBFlt*) SafeMalloc( numBranches * nPairs * sizeof(MrBFlt));
            if (D1LP[i] == NULL)
                return (ERROR);
            }

        /*  init arrays for hessian and jacobian */
        H =    (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        J =    (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        Hinv = (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        HiJ =  (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        for (i=0; i<numBranches; i++) 
            { 
            H[i]=(MrBFlt*) SafeMalloc(   numBranches * sizeof(MrBFlt));
            J[i]=(MrBFlt*) SafeMalloc(   numBranches * sizeof(MrBFlt));
            Hinv[i]=(MrBFlt*) SafeMalloc(numBranches * sizeof(MrBFlt));
            HiJ[i]=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
            }

        eigvals=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
        if (eigvals == NULL)
                return (ERROR);
        eigvalsc=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
        if (eigvalsc == NULL)
                return (ERROR);

        /* alloc helper vector for storing partition pairs   */
        tempPartitionPair=(int*)SafeMalloc(numLocalTaxa * sizeof(int));

        /*  *
         *  Done initializing
         *  */

        /* first calculate the n_ii for the division */
        splitI=0;
        for (k=0; k<numLocalTaxa-1; k++) 
            {
            for (l=k+1; l<numLocalTaxa; l++)
                {
                pI=pairIdx(k,l,numLocalTaxa); /*  just get the single pair idx */
                overallPwIdx=countIndex[nSplits][pI];
                for (c=0; c<numChar; c++)
                    {
                    if (partitionId[c][partitionNum] != d+1) /* only count within partition */
                        continue; 

                    cI=countIndex[splitI][pI];

                    c1=toIdx(matrix[pos(k,c,numChar)]);
                    c2=toIdx(matrix[pos(l,c,numChar)]);

                    if (c1 < 0 || c2 < 0) 
                        continue;
                    dI=dIdx(c1,c2,nStates);

                    counts[cI+dI] += 1;
                    counts[overallPwIdx+dI] += 1;
                    splitI+=1;
                    splitI=splitI%nSplits;

                    if (cI+dI < 0 || overallPwIdx+dI < 0)
                            MrBayesPrint("possible oob? %d %d %d", cI, overallPwIdx, dI);
                    }
                }
            }

        /*  get the counts per data split  */
        for (k=0; k<nPairs; k++) 
            {   
            for (splitI=0; splitI<nSplits; splitI++) 
                {
                cI=countIndex[splitI][k];
                nI=niiIndex[splitI][k];
                nIOverall=niiIndex[nSplits][k];
                for (i=0;i<nStates;i++) 
                    {
                    for (j=0;j<nStates;j++) 
                        {
                        dI = dIdx(i,j,nStates);
                        count=counts[cI+dI];
                        if (i==j) 
                            {
                            n11[nI]        += count;
                            n11[nIOverall] += count; 
                            }
                        else 
                            {
                            n10[nI]        += count;
                            n10[nIOverall] += count;
                            }
                        }
                    }
                }
            }

        /*  now just calculate the pw dists */
        
        if (m->shape != NULL)
            al=*GetParamVals(m->shape,chain,state[chain]);
        else 
            al=0.0;

        for (k=0; k<nPairs; k++) 
            {   
            nI=niiIndex[nSplits][k];
            MrBFlt prop =  (n10[nI] * 1.0) / (n10[nI] + n11[nI]);
            if (n10[nI] == 0) 
                pwDists[k] = 0.0;
            else if (n11[nI] == 0)
                pwDists[k] = TIME_MAX; /* equilibrium...  */
            else if (prop > 0.75) 
                pwDists[k] = TIME_MAX;
            else 
                {
                if (m->shape != NULL)
                    pwDists[k] = al * (3.0/4) * (pow(1-(4 * prop/3), -(1.0/al)) - 1.0);
                else 
                    pwDists[k] = -1.0 * (3.0/4) * log(1.0 - (4.0/3.0) * prop);
                }
                if (isnan(pwDists[k])) 
                    MrBayesPrint("nan dist\n");
                //MrBayesPrint(" %f \n", pwDists[k]);

            } 

        /* get base rate */
        baseRate = GetRate (d, chain);
    
        /* compensate for invariable sites if appropriate */
        if (m->pInvar != NULL)
            baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
       
        /* get category rates */
        theRate = 1.0;
        if (m->shape != NULL)
            catRate = GetParamSubVals (m->shape, chain, state[chain]);
        else if (m->mixtureRates != NULL)
            catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
        else
            catRate = &theRate;

        /*  set up pair/branch indicator matrix:  */
        nLongsNeeded=((numLocalTaxa-1)/nBitsInALong)+1;

        PairBranch = SafeMalloc( numBranches * sizeof(int*)) ;
        if (!PairBranch) 
            return(ERROR);

        for (i=0; i<numBranches; i++) 
            {
            PairBranch[i] = SafeMalloc(nPairs * sizeof(int));
            if(!PairBranch[i])
                return(ERROR);
            }

        l1=0;
        l2=0;
        TreeNode *p;

        // Make sure we have bitfields allocated and set
        if (tree->bitsets == NULL)
            {
            AllocateTreePartitions(tree);
            freeBitsets = YES;
            }
        else
            {
            ResetTreePartitions(tree);   // just in case
            freeBitsets = NO;
            }

        for (i=index=0;i<tree->nNodes;i++)
            {
            p=&(tree->nodes[i]);
            if (AreDoublesEqual(p->length,0.0,ETA)) continue;

            for (j=0;j<numLocalTaxa;j++) /*  reset helper array */
                tempPartitionPair[j]=1;

            for (l1=FirstTaxonInPartition(p->partition, nLongsNeeded); 
                 l1<numLocalTaxa; 
                 l1=NextTaxonInPartition(l1, p->partition, nLongsNeeded))
                 tempPartitionPair[l1]=0; /*  now temp array has 1s for taxa not in partition */

            /*  now loop again and fill in 1s for pairs with a taxa in this partition and one not in partition */
            for (l2=FirstTaxonInPartition(p->partition, nLongsNeeded); 
                 l2<numLocalTaxa; 
                 l2=NextTaxonInPartition(l2, p->partition, nLongsNeeded))
                { 
                for (j=0;j<numLocalTaxa;j++)
                    {
                    if (j == l2) continue;
                    if (tempPartitionPair[j] == 1) 
                        {
                        k=pairIdx(l2,j,numLocalTaxa);
                        PairBranch[index][k]=1; /*  both taxa below node, so node not in pair path  */
                        }
                    }
                }
            index++; 
            }

        /*  calculate derivatives needed for J/H */
        for (k=0; k<nPairs; k++) 
            {
            
            nRates=m->numRateCats; 
            if (pwDists[k] == 0.0) 
                {
                p_10[k] = 0.25;
                p_11[k] = 0.75;
                p1_10[k] = 0.0;
                p1_11[k] = 0.0;
                }
            else if (pwDists[k] >= TIME_MAX)
                {
                p_10[k] = 0.0;
                p_11[k] = 1.0;
                p1_10[k] = 0.0;
                p1_11[k] = 0.0;
                }
            else 
                {
                for (l=0; l<nRates; l++)
                    {
                    rc =  baseRate * catRate[l];
                    eterm = exp(-(4.0/3)*rc*pwDists[k]);
                    p_10[k]  +=  (1.0/nRates) * ((3.0/4) - (3.0/4) * eterm);
                    p_11[k]  +=  (1.0/nRates) * ((1.0/4) + (3.0/4) * eterm);
                    p1_10[k] += ( 1.0 / nRates) * rc * eterm;
                    p1_11[k] += (-1.0 / nRates) * rc * eterm;
                    //p2_10[k] += (-4.0 / (3.0*nRates)) * rc*rc * eterm;
                    //p2_11[k] += ( 4.0 / (3.0*nRates)) * rc*rc * eterm;
                    }
                }
            }

        //for (k=0;k<nPairs;k++)
        //    {
        //    MrBayesPrint(" -- pair %d -- \n", k);
        //    MrBayesPrint(" %f \n", pwDists[k]);
        //    MrBayesPrint(" %f \n", p_10[k]);
        //    MrBayesPrint(" %f \n", p_11[k]);
        //    MrBayesPrint(" %f \n", p1_10[k]);
        //    MrBayesPrint(" %f \n", p1_11[k]);
        //    }

        /*  Now compute first derivs of composite ll derivs  */
        for (i=index=0; i<numBranches; i++)
            {
            for (k=0; k<nPairs; k++)
                {

                if (p_11[k] < ETA) 
                    t1=0.0;
                else 
                    t1 = (p1_11[k] / p_11[k]);

                if (p_10[k] < ETA)
                    t2=0.0;
                else 
                    t2 = (p1_10[k] / p_10[k]);

                for (z=0; z<nSplits; z++)
                    {
                    nidx=niiIndex[z][k];


                    if (PairBranch[i][k]==1) 
                        {
                        /*  fill in derivative arrays */
                        //MrBayesPrint("%d %d %d  %d, %f \n", i, k, z, n11[nidx], t1);
                        //MrBayesPrint("%d %d %d  %d, %f \n", i, k, z, n10[nidx], t2);
                        //MrBayesPrint("\n");

                        if (isnan((n11[nidx] * t1 + n10[nidx] * t2)))
                            MrBayesPrint("NaN\n");

                        D1L[z][i]     += (n11[nidx] * t1 + n10[nidx] * t2);
                        D1LP[z][index] = (n11[nidx] * t1 + n10[nidx] * t2);

                        }
                    }
                index++;
                }
            }

        // Debug output
        // MrBayesPrint("D1L\n");
        // for (i=0;i<4;i++) {
        //     for (j=0;j<4;j++) {
        //         MrBayesPrint(" %f ", D1L[i][j]);
        //     }
        //     MrBayesPrint("\n");
        // }
        //
        // MrBayesPrint("D1LP\n");
        // for (i=0;i<4;i++) {
        //     for (j=0;j<4;j++) {
        //         MrBayesPrint(" %f ", D1LP[i][j]);
        //     }
        //     MrBayesPrint("\n");
        // }
           
        // MrBayesPrint(" %d ", n11[niiIndex[0][0]]);
        // MrBayesPrint(" %d \n", n10[niiIndex[0][0]]);
        //
        // MrBayesPrint(" %d ", n11[niiIndex[4][0]]);
        // MrBayesPrint(" %d \n", n10[niiIndex[4][0]]);

        /*  fill in J and H  */
        for (i=0; i<numBranches; i++)
            {
            for (j=i; j<numBranches; j++)
                {
                         
                H[i][j]=0.0;
                J[i][j]=0.0;

                for (z=0; z<nSplits; z++)
                    {
                    //nidx=niiIndex[d][k];
                    J[i][j] += (1.0/nSplits) * D1L[z][i] * D1L[z][j] ;
                    for (k=0; k<nPairs; k++ )
                        {
                        if (PairBranch[i][k]==1 && PairBranch[j][k]==1)
                            {
                            if (nPairs*i+k >= nPairs * numBranches || nPairs*j+k >= nPairs * numBranches)
                                {
                                    MrBayesPrint("possible index oob: %d , %d, %d ", nPairs*j+k, nPairs*i+k, nPairs * numBranches );
                                }
                            H[i][j] += (1.0/nSplits) * D1LP[z][nPairs*i+k] * D1LP[z][nPairs*j+k] ;
                            }
                        }
                    }

                if (i != j) 
                    {
                    J[j][i]=J[i][j];
                    H[j][i]=H[i][j];
                    }

                }
            }

        /*  compute  H^-1 * J and the eigenvalues:  */
        dw= (MrBFlt *)SafeMalloc((size_t)numBranches*(sizeof(MrBFlt)));
        iw= (int *)SafeMalloc((size_t)numBranches*(sizeof(int)));

 
        InvertMatrix(numBranches, H, dw,iw, Hinv);
        MultiplyMatrices(numBranches, Hinv, J, HiJ);


        isComplex=GetEigens(numBranches,HiJ,eigvals,eigvalsc,V,Vinv,Vc,Vcinv);

        //MrBayesPrint("isComplex: %d \n", isComplex);

        eigsum=0.0;
        eigsum2=0.0 ;
        em=0.0;
        v=0.0;

        for (i=0; i<numBranches; i++) {
            // MrBayesPrint("Eigen %d = %f \n", i, eigvals[i]);
            eigsum += fabs(eigvals[i]);
            eigsum2 += eigvals[i] * eigvals[i];
        }

        em = eigsum/(1.0*(numBranches)); 
        v = (eigsum) / eigsum2;

        if (m->usePwWeights == 1)
            pwWeight[d]=(1.0) / em;
        else if (m->usePwWeights == 2)  
            pwWeight[d]=v;

        MrBayesPrint("%s PW Weight (division %d): %f \n", spacer, d, pwWeight[d]);

        /*  free allocations   */
        /*  helper matrices */

        //MrBayesPrint("%s 1st chunk \n", spacer);        
        FreeSquareDoubleMatrix(V);
        FreeSquareDoubleMatrix(Vinv);
        FreeSquareComplexMatrix(Vc);
        FreeSquareComplexMatrix(Vcinv);
        free(tempPartitionPair);

        free(dw);
        free(iw);

        /*  counts */
        //MrBayesPrint("%s 2nd chunk \n", spacer);        
        free(counts);
        //MrBayesPrint("%s 2a \n", spacer);        
        for (i=0; i<(nSplits+1); i++)
            free(countIndex[i]);
        //MrBayesPrint("%s 2b \n", spacer);        
        free(countIndex);
        for (i=0; i<(nSplits+1); i++)
            free(niiIndex[i]);
        //MrBayesPrint("%s 2c \n", spacer);        
        free(niiIndex);
        //MrBayesPrint("%s 2d \n", spacer);        
        free(n10);
        //MrBayesPrint("%s 2e \n", spacer);        
        free(n11);

        /*  dists and probabilities */
        //MrBayesPrint("%s 1st free \n", spacer);        
        free(pwDists);
        //MrBayesPrint("%s 2nd free \n", spacer);
        free(p_10 );
        //MrBayesPrint("%s 3nd free \n", spacer);
        free(p_11 );
        //MrBayesPrint("%s 4th free \n", spacer);
        free(p1_10); 
        //MrBayesPrint("%s 5th free \n", spacer);
        free(p1_11);  
                      
        /*  hessian and jacobian */
        for (i=0; i<numBranches; i++) 
            { 
            free(H[i]);
            free(J[i]);
            free(HiJ[i]);
            free(Hinv[i]);
            }
        free(H);
        free(J);
        free(Hinv);
        free(HiJ);

        for (i=0; i<nSplits; i++) 
            { 
            free(D1L[i]);
            free(D1LP[i]);
            }
        free(D1L);
        free(D1LP);
        free(eigvals);
        free(eigvalsc);

        //MrBayesPrint("%s 6th free \n", spacer);
        for (i=0; i<numBranches; i++)
            free(PairBranch[i]);
        free(PairBranch);

        } /* end loop over numCurrentDivisions */

    return(0);
}

MrBFlt EstPwDist_GTR(ModelInfo *m, int chain, int* counts, int countIdx, MrBFlt al)
{
    int i,j,pid;
    MrBFlt tau;
    MrBFlt *bs;

    MrBFlt **V;
    MrBFlt **Vinv;      
    MrBFlt **LaLog;     
                       
    MrBComplex **Vc;    
    MrBComplex **Vcinv; 
    MrBFlt la[4]; 
    MrBFlt laC[4];
                       
    MrBFlt **F;            
    MrBFlt **Q;            
    MrBFlt **Temp;         

    // set up matrices 
    V     = AllocateSquareDoubleMatrix(4);
    Vinv  = AllocateSquareDoubleMatrix(4);
    LaLog = AllocateSquareDoubleMatrix(4);

    Vc    = AllocateSquareComplexMatrix(4); 
    Vcinv = AllocateSquareComplexMatrix(4);

    F = AllocateSquareDoubleMatrix(4);
    Q = AllocateSquareDoubleMatrix(4);
    Temp = AllocateSquareDoubleMatrix(4);
    
    tau=0.0;

    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    // calc denom:
    int tot=0;
    for (i=0;i<4;i++) 
        for (j=0;j<4;j++)
            tot += counts[countIdx + dIdx(i,j,4)];

    // set up matrix of empirical transitions 
    //   probabilities, 'Q' (will be modified in place)  
    for (i=0;i<4;i++) {
        for (j=0;j<4;j++){
            F[i][j] = 1.0 * counts[countIdx + dIdx(i,j,4)] / (tot); 
        }
    }

    // Symmetrize F:
    for (i=0;i<4;i++) {
        for (j=i;j<4;j++){
            F[i][j] = (F[i][j]+F[j][i])/2.0;
            if (i != j)
                F[j][i] = F[i][j];
        }
    }

    // need to multiply rows by rowsums of F 
    MrBFlt rowSums[4] = {0.0};
    for (i=0;i<4;i++) 
        for (j=0;j<4;j++)
            rowSums[i] += F[i][j];
    for (i=0;i<4;i++) 
        for (j=0;j<4;j++)
            F[i][j] = F[i][j] / rowSums[i]; 

    int isComplex=GetEigens(4,F,la,laC,V,Vinv,Vc,Vcinv);
    (void)isComplex;
   
    // diagonal matrix of log^{lambda_i}, lambdas are eigvals of Q
    for (i=0;i<4;i++) {
        if (al > 0.0) 
            LaLog[i][i]= al * (1- pow(la[i], -(1.0/al)));
        else 
            LaLog[i][i]=log(la[i]);

        for (j=0;j<i;j++) {
            LaLog[i][j]=0.0;
            LaLog[j][i]=0.0;
        }
    }

    // calculate the matrix log: log(e^Qt) = V log[La] V^-1) 
    MultiplyMatrices(4,V,LaLog,Temp); 
    MultiplyMatrices(4,Temp,Vinv,Q); // Q is ptr to resulting matrix
    
    // set diagonal so rowsums are 0, mult by inverse Dpi mat:
    double rowsum;
    for (i=0;i<4;i++) {
        rowsum=0.0;
        for (j=0;j<4;j++) {
            if (j!=i) rowsum+=Q[i][j];
        }
        Q[i][i]=-1.0*rowsum;
    }
    
    for (i=0;i<4;i++)
        tau += -1.0 * Q[i][i] * rowSums[i];

    // set up matrices 
    FreeSquareDoubleMatrix(V);
    FreeSquareDoubleMatrix(Vinv);
    FreeSquareDoubleMatrix(LaLog);

    FreeSquareComplexMatrix(Vc); 
    FreeSquareComplexMatrix(Vcinv);

    FreeSquareDoubleMatrix(F);
    FreeSquareDoubleMatrix(Q);
    FreeSquareDoubleMatrix(Temp);
  
    return(tau);
}


int TranProbMatrix_GTR(ModelInfo *m, int chain, double dist, double al, double *transProbs)
{
    int i,j,k;

    MrBFlt **V;
    MrBFlt **Vinv;      
    MrBFlt **LaLog;     
                       
    MrBComplex **Vc;    
    MrBComplex **Vcinv; 
    MrBFlt la[4]; 
    MrBFlt laC[4];
                       
    MrBFlt **F;            
    MrBFlt **Q;            
    MrBFlt **Temp;         
    MrBFlt **TransProbTemp;         

    MrBFlt *bs;
    MrBFlt *rateValues;
    MrBFlt *catRate;

    MrBFlt **LaExp;   
    MrBFlt **TempMat;   
    MrBFlt **Qtausr;   

    // set up matrices 
    V     = AllocateSquareDoubleMatrix(4);
    Vinv  = AllocateSquareDoubleMatrix(4);
    LaLog = AllocateSquareDoubleMatrix(4);

    Vc    = AllocateSquareComplexMatrix(4); 
    Vcinv = AllocateSquareComplexMatrix(4);

    F = AllocateSquareDoubleMatrix(4);
    Q = AllocateSquareDoubleMatrix(4);
    Temp = AllocateSquareDoubleMatrix(4);
    TransProbTemp = AllocateSquareDoubleMatrix(4);
    
    // compute gtr transition probabilities
    // set up matrices for taking the matrix exponent
    LaExp= AllocateSquareDoubleMatrix(4);
    TempMat= AllocateSquareDoubleMatrix(4);
    Qtausr= AllocateSquareDoubleMatrix(4);

    // set up q matrix
    rateValues = GetParamVals(m->revMat, chain, state[chain]);
    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    if (al > 0.0)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);

    /* reset input matrix*/
    for (i=0;i<4;i++)
       for (j=0;j<4;j++)
           transProbs[dIdx(i,j,4)] = 0.0;

    /* set diagonal of Q matrix to 0 */
    for (i=0; i<4; i++)
        Q[i][i] = 0.0;
  
    /* initialize Q matrix */
    MrBFlt scaler, mult;
    scaler = 0.0;
    for (i=0; i<4; i++)
        {
        for (j=i+1; j<4; j++)
            {
            if (i == 0 && j == 1)
                mult = rateValues[0];
            else if (i == 0 && j == 2)
                mult = rateValues[1];
            else if (i == 0 && j == 3)
                mult = rateValues[2];
            else if (i == 1 && j == 2)
                mult = rateValues[3];
            else if (i == 1 && j == 3)
                mult = rateValues[4];
            else if (i == 2 && j == 3)
                mult = rateValues[5];
            Q[i][i] -= (Q[i][j] = bs[j] * mult);
            Q[j][j] -= (Q[j][i] = bs[i] * mult);
            scaler += bs[i] * Q[i][j];
            scaler += bs[j] * Q[j][i];
            }
        }
       
    /* rescale Q matrix */
    scaler = 1.0 / scaler;
    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
            Q[i][j] *= scaler;
    
    // now compute the transition probability matrix:
    for (k=0; k<m->numRateCats; k++)
        {
        // probability transition matrix for site rate i:
        if (al > 0.0) {
            MultiplyMatrixByScalar(4, Q, dist * catRate[k], Qtausr);  
        }
        else 
            MultiplyMatrixByScalar(4, Q, dist, Qtausr);  

        int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
           
        // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
        for (i=0;i<4;i++) {
            LaExp[i][i] = exp(la[i]);
            for (j=0;j<i;j++) {
                LaExp[i][j]=0.0;
                LaExp[j][i]=0.0;
            }
        }
           
        // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
        MultiplyMatrices(4,V,LaExp,TempMat); 
        MultiplyMatrices(4,TempMat,Vinv,TransProbTemp);

        for (i=0;i<4;i++)
            for (j=0;j<4;j++)
                transProbs[dIdx(i,j,4)] += bs[i] * (1.0/m->numRateCats) * TransProbTemp[i][j];

        //.  turn into doublet probs */
        //for (i=0;i<4;i++)
        //    for (j=0;j<4;j++)
        //        transProbs[dIdx(i,j,4)] *= bs[i];

    }

    // free matrices 
    FreeSquareDoubleMatrix(V);
    FreeSquareDoubleMatrix(Vinv);
    FreeSquareDoubleMatrix(LaLog);

    FreeSquareComplexMatrix(Vc); 
    FreeSquareComplexMatrix(Vcinv);

    FreeSquareDoubleMatrix(F);
    FreeSquareDoubleMatrix(Q);
    FreeSquareDoubleMatrix(Temp);
    FreeSquareDoubleMatrix(TransProbTemp);
    
    // compute gtr transition probabilities
    // set up matrices for taking the matrix exponent
    FreeSquareDoubleMatrix(LaExp);
    FreeSquareDoubleMatrix(TempMat);
    FreeSquareDoubleMatrix(Qtausr);

    return(1);
}

int CalcPairwiseWeights_GTR (int chain) {

    /*  setup the arays for pairwise counts,
     *  and populate the counts.
     * 
     *  this function will supercede the 'PairwiseCounts'
     *  function, since that just uses a global variable and 
     *  we want to be able to apply pw likelihood within partitions.
     *
     *  We'll also make and populate the counts for the data splits
     *  for estimating jacobians for weighting the pw likelihood. 
     *  */

    ModelInfo* m;

    int i,j,k,l,c,c1,c2,d,z;
    int pI,cI,dI,splitI;
    int nI;
    int *counts;
    int **countIndex;
    int overallPwIdx;
    int index, indexStep;
    int nSplits, nPairs, nStates;
    MrBFlt *pwDists;
    MrBFlt baseRate;
    MrBFlt *catRate;
    MrBFlt theRate;
    //MrBFlt rc;
    MrBFlt **tp, **tp1; 
    MrBFlt *tptemp;
    MrBFlt t1;

    int    nidx;
    MrBFlt **J, **H, **Hinv;
    MrBFlt **HiJ, *eigvals, *eigvalsc;
    int **niiIndex ;
    int numBranches;
    Tree *tree;
    int freeBitsets;
    int *tempPartitionPair;

    //MrBFlt h=MRBFLT_MIN;
    //MrBFlt h=ETA;
    MrBFlt h=1E-10;
    MrBFlt *bs;

    /*  first set up worker matrices for eigen computation */
    // set up matrices 
    MrBFlt **V, **Vinv;
    MrBComplex **Vc, **Vcinv; 

    // loop over the partitions:
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];

        bs = GetParamSubVals(m->stateFreq, chain, state[chain]);
        tree = GetTree(m->brlens, chain, state[chain]);
        if (m->usePwWeights==3)
            {
            pwWeight[d] = 2.0 / (numLocalTaxa * (numLocalTaxa - 1));
            MrBayesPrint("%s pwWeight: %f \n", spacer, m->pwWeight);
            continue;
            }
                
        nStates = m->numModelStates;
        nSplits = m->numDataSplits;
        nPairs = m->numPairs;

        if (!tree->isRooted)
            numBranches=numLocalTaxa*2 - 3;
        else 
            numBranches=numLocalTaxa*2 - 2;

        overallPwIdx = m->numDataSplits;
        MrBayesPrint("%s Calculating PW Weight using method %f \n",  spacer, m->usePwWeights );

        /*  * 
         *  Initialize necessary arrays:
         *  */
        /* initialize counts & index array */

        V     = AllocateSquareDoubleMatrix(numBranches);
        Vinv  = AllocateSquareDoubleMatrix(numBranches);
             
        Vc    = AllocateSquareComplexMatrix(numBranches); 
        Vcinv = AllocateSquareComplexMatrix(numBranches);


        int countLen = (nSplits+1) * nPairs * nStates * nStates;
        // MrBayesPrint("%s count array length: %d \n", spacer, countLen);

        counts = (int*) SafeMalloc( (nSplits+1) * nPairs * nStates * nStates * sizeof(int));
        if (counts == NULL)
            return(ERROR);

        // MrBayesPrint("%s Done alloc counts \n", spacer);
        countIndex=(int**) SafeMalloc( (nSplits+1) * sizeof(int*));
        if (countIndex == NULL)
            return(ERROR);        

        for (i=0; i<(nSplits+1); i++)
            {
            countIndex[i]=(int*) SafeMalloc( nPairs * sizeof(int));
            if (countIndex[i] == NULL)
                return(ERROR);        
            }


        index=0;
        indexStep=nStates*nStates;
        for (i=0; i<(nSplits+1); i++)
            {
            for (j=0; j<nPairs; j++)
                {
                countIndex[i][j]=index;
                index+=indexStep;
                }
            }

        /*  init array for pw distances */
        pwDists = (MrBFlt*) SafeMalloc( nPairs * sizeof(MrBFlt));
        if (pwDists == NULL) 
            return (ERROR);

        /* init arrays for derivatives  */
        /*  TODO: each deriv will be 4x4 matrix (Q * exp(tau * Q)) */
        /*  p is transition probabilities, p1 is 1st derive of transition probs */
        tp = (MrBFlt**) SafeMalloc( nPairs * sizeof(MrBFlt*));
        if (tp == NULL)
            return(ERROR);

        tp1 = (MrBFlt**) SafeMalloc( nPairs * sizeof(MrBFlt*));
        if (tp1 == NULL)
            return(ERROR);

        for (k=0;k<nPairs;k++)
            {
            tp[k] = (MrBFlt*) SafeMalloc(16 * sizeof(MrBFlt));
            if (tp[k] == NULL)
                return(ERROR);

            tp1[k] = (MrBFlt*) SafeMalloc(16 * sizeof(MrBFlt));
            if (tp1[k] == NULL)
                return(ERROR);
            }

        tptemp = (MrBFlt*) SafeMalloc( 16 * sizeof(MrBFlt));
        if (tptemp == NULL)
            return(ERROR);


        /*  init arrays for first and second derivs */
        MrBFlt  **D1L, **D1LP;
        D1L =  (MrBFlt**) SafeMalloc( nSplits * sizeof(MrBFlt*));
        if (!D1L)
            return (ERROR);

        for (i=0; i<nSplits; i++) 
            { 
            D1L[i]=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
            if (!D1L[i])
                return (ERROR);
            }

        D1LP =  (MrBFlt**) SafeMalloc( nSplits * sizeof(MrBFlt*));
        if (!D1LP)
            return (ERROR);

        for (i=0; i<nSplits; i++) 
            { 
            D1LP[i]=(MrBFlt*) SafeMalloc( numBranches * nPairs * sizeof(MrBFlt));
            if (D1LP[i] == NULL)
                return (ERROR);
            }

        /*  init arrays for hessian and jacobian */
        H =    (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        J =    (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        Hinv = (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        HiJ =  (MrBFlt**) SafeMalloc( numBranches * sizeof(MrBFlt*));
        for (i=0; i<numBranches; i++) 
            { 
            H[i]=(MrBFlt*) SafeMalloc(   numBranches * sizeof(MrBFlt));
            J[i]=(MrBFlt*) SafeMalloc(   numBranches * sizeof(MrBFlt));
            Hinv[i]=(MrBFlt*) SafeMalloc(numBranches * sizeof(MrBFlt));
            HiJ[i]=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
            }

        eigvals=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
        if (eigvals == NULL)
                return (ERROR);
        eigvalsc=(MrBFlt*) SafeMalloc( numBranches * sizeof(MrBFlt));
        if (eigvalsc == NULL)
                return (ERROR);

        /* alloc helper vector for storing partition pairs   */
        tempPartitionPair=(int*)SafeMalloc(numLocalTaxa * sizeof(int));

        /*  *
         *  Done initializing
         *  */

        /* first calculate the nucleotide pair counts for the division */
        splitI=0;
        for (k=0; k<numLocalTaxa-1; k++) 
            {
            for (l=k+1; l<numLocalTaxa; l++)
                {
                pI=pairIdx(k,l,numLocalTaxa); /*  just get the single pair idx */
                overallPwIdx=countIndex[nSplits][pI];
                for (c=0; c<numChar; c++)
                    {
                    if (partitionId[c][partitionNum] != d+1) /* only count within partition */
                        continue; 

                    //splitI=(c*nSplits)/numChar;
                    cI=countIndex[splitI][pI];

                    c1=toIdx(matrix[pos(k,c,numChar)]);
                    c2=toIdx(matrix[pos(l,c,numChar)]);

                    if (c1 < 0 || c2 < 0) 
                        continue;

                    dI=dIdx(c1,c2,nStates);
                    counts[cI + dI] += 1;
                    counts[overallPwIdx + dI] += 1;
                    if (c > splitI * (numChar / nSplits)) splitI += 1;
                    // splitI+=1;
                    splitI = splitI % nSplits;

                    if (cI+dI < 0 || overallPwIdx+dI < 0)
                        MrBayesPrint("possible oob? %d %d %d", cI, overallPwIdx, dI);
                    }
                }
            }

        /*  now just calculate the pw dists */
        MrBFlt al;
        if (m->shape != NULL)
            al=*GetParamVals(m->shape,chain,state[chain]);
        else 
            al=0.0;


        /*  can use EstPwDist_GTR function */
        for (k=0; k<nPairs; k++) 
            {   
            nI=countIndex[nSplits][k];
            pwDists[k] = EstPwDist_GTR(m, chain, counts, nI, al);
            if (isnan(pwDists[k])) 
                MrBayesPrint("nan dist\n");
            }

        /* get base rate */
        baseRate = GetRate (d, chain);
    
        /* compensate for invariable sites if appropriate */
        if (m->pInvar != NULL)
            baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
       
        /* get category rates */
        theRate = 1.0;
        if (m->shape != NULL)
            catRate = GetParamSubVals (m->shape, chain, state[chain]);
        else if (m->mixtureRates != NULL)
            catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
        else
            catRate = &theRate;

        /*  set up pair/branch indicator matrix:  */
        int nLongsNeeded=((numLocalTaxa-1)/nBitsInALong)+1;

        int **PairBranch;
        PairBranch = SafeMalloc( numBranches * sizeof(int*)) ;
        if (!PairBranch) 
            return(ERROR);

        for (i=0; i<numBranches; i++) 
            {
            PairBranch[i] = SafeMalloc(nPairs * sizeof(int));
            if(!PairBranch[i])
                return(ERROR);
            }

        int l1;
        int l2;
        TreeNode *p;

        // Make sure we have bitfields allocated and set
        if (tree->bitsets == NULL)
            {
            AllocateTreePartitions(tree);
            freeBitsets = YES;
            }
        else
            {
            ResetTreePartitions(tree);   // just in case
            freeBitsets = NO;
            }

        for (i=index=0;i<tree->nNodes;i++) // set for each node in tree
            {
            p=&(tree->nodes[i]);
            if (AreDoublesEqual(p->length,0.0,ETA)) continue;

            for (j=0;j<numLocalTaxa;j++) /*  reset helper array */
                tempPartitionPair[j]=1;

            for (l1=FirstTaxonInPartition(p->partition, nLongsNeeded); 
                 l1<numLocalTaxa; 
                 l1=NextTaxonInPartition(l1, p->partition, nLongsNeeded))
                 tempPartitionPair[l1]=0; /*  now temp array has 1s for taxa not in partition */

            /*  now loop again and fill in 1s for pairs with a taxa in this partition and one not in partition */
            for (l2=FirstTaxonInPartition(p->partition, nLongsNeeded); 
                 l2<numLocalTaxa; 
                 l2=NextTaxonInPartition(l2, p->partition, nLongsNeeded))
                { 
                for (j=0;j<numLocalTaxa;j++)
                    {
                    if (j == l2) continue;
                    if (tempPartitionPair[j] == 1) 
                        {
                        k=pairIdx(l2,j,numLocalTaxa);
                        PairBranch[index][k]=1;     /*  l1 in and l2 not in partition, so path btw l1,lw split by node  */
                        }
                    }
                }
            index++; 
            }

        /*  calculate derivatives needed for J/H */
        /*  This needs to change for GTR model -- compute transition probability matrix 
         *  derivatives numerically...  */

        for (k=0; k<nPairs; k++) 
            {
            if (pwDists[k] == 0.0) 
                {
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        z=dIdx(i,j,4);
                        if (i == j)
                            tp[k][z] = 1;
                        else 
                            tp[k][z] = 0;
                        tp1[k][z] = 0;
                        }
                    }
                }
            else if (pwDists[k] >= TIME_MAX)
                {
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        z=dIdx(i,j,4);
                        tp[k][z] = bs[i];
                        tp1[k][z] = 0;
                        }
                    }
                }
            else 
                {
                /*  compute transition probabilities and numerical deriv  */
                TranProbMatrix_GTR(m, chain, pwDists[k], al, tp[k]);

                // numerical derivative
                TranProbMatrix_GTR(m, chain, pwDists[k]+h, al, tptemp);
                for (i=0;i<4;i++)
                    for (j=0;j<4;j++)
                        tp1[k][dIdx(i,j,4)] = (tptemp[dIdx(i,j,4)] - tp[k][dIdx(i,j,4)]) / h  ;
                }
            }

        //for (k=0;k<nPairs;k++)
        //    {
        //    MrBayesPrint("--Pair %d --\n",k  );
        //    MrBayesPrint("%f \n", pwDists[k]);

        //    for (i=0;i<16;i++)
        //        MrBayesPrint("%f ", tp[k][i]);
        //    MrBayesPrint("\n");

        //    for (i=0;i<16;i++)
        //        MrBayesPrint("%f ", tp1[k][i]);
        //    MrBayesPrint("\n");
        //    MrBayesPrint("\n");
        //    }

        /*  Now compute first derivs of composite ll   */
        for (i=index=0; i<numBranches; i++)
            {
            for (k=0; k<nPairs; k++)
                {
                if (PairBranch[i][k] == 1) 
                    {
                    for (j=0;j<16;j++)
                        {
                        if (tp[k][j] < ETA) 
                            continue;

                        t1 = (tp1[k][j] / (tp[k][j]));
                        for (z=0; z<nSplits; z++)
                            {
                            nidx=countIndex[z][k];

                            /*  fill in derivative arrays */
                            //if (isnan(counts[nidx+j] * t1)) {
                            //    continue;
                            //}

                            D1L[z][i]      += (counts[nidx+j] * t1);
                            D1LP[z][index] += (counts[nidx+j] * t1);

                            }
                        }
                    }
                index++;
                }
            }

        //MrBayesPrint("D1L\n");
        //for (i=0;i<4;i++) {
        //    for (j=0;j<4;j++) {
        //        MrBayesPrint(" %f ", D1L[i][j]);
        //    }
        //    MrBayesPrint("\n");
        //}

        //MrBayesPrint("D1LP\n");
        //for (i=0;i<4;i++) {
        //    for (j=0;j<4;j++) {
        //        MrBayesPrint(" %f ", D1LP[i][j]);
        //    }
        //    MrBayesPrint("\n");
        //} 
        /*  fill in J and H  */

         for (i=0; i<numBranches; i++)
             {
             for (j=i; j<numBranches; j++)
                 {

               H[i][j]=0.0;
               J[i][j]=0.0;

               for (z=0; z<nSplits; z++)
                   {
                   J[i][j] += (1.0/nSplits) * D1L[z][i] * D1L[z][j] ;
                   for (k=0; k<nPairs; k++ )
                       {
                       if (PairBranch[i][k]==1 && PairBranch[j][k]==1)
                           {
                           if (nPairs*i+k >= nPairs * numBranches || nPairs*j+k >= nPairs * numBranches)
                                 MrBayesPrint("possible index oob: %d , %d, %d ", nPairs*j+k, nPairs*i+k, nPairs * numBranches );
                           H[i][j] += (1.0/nSplits) * D1LP[z][nPairs*i + k] * D1LP[z][nPairs*j + k] ;
                           }
                       }
                   }

                if (i != j) 
                    {
                    J[j][i]=J[i][j];
                    H[j][i]=H[i][j];
                    }

                }
            }

        /*  compute  H^-1 * J and the eigenvalues:  */
        MrBFlt* dw= (MrBFlt *)SafeMalloc((size_t)numBranches*(sizeof(MrBFlt)));
        int*    iw= (int *)SafeMalloc((size_t)numBranches*(sizeof(int)));

        InvertMatrix(numBranches, H, dw,iw, Hinv);
        MultiplyMatrices(numBranches, Hinv, J, HiJ);

        int isComplex=GetEigens(numBranches,HiJ,eigvals,eigvalsc,V,Vinv,Vc,Vcinv);
        //for (k=0; k<numBranches; k++)
        //    {
        //    MrBayesPrint("HiJ[%d][.] =",k);
        //    for (l=0; l<numBranches; l++)
        //        {
        //        MrBayesPrint("  % .3f", HiJ[k][l]);
        //        }
        //        MrBayesPrint("\n");
        //    }

        MrBFlt eigsum=0.0 ;
        MrBFlt eigsum2=0.0 ;
        MrBFlt em=0.0;
        MrBFlt v=0.0;

        for (i=0; i<numBranches; i++) {
            //MrBayesPrint("Eigen %d = %f \n", i, eigvals[i]);
            eigsum += fabs(eigvals[i]) ;
            eigsum2 += eigvals[i] * eigvals[i];
        }

        em = eigsum/(1.0*(numBranches)); 
        v = eigsum / eigsum2;

        if (m->usePwWeights == 1)
            pwWeight[d]=(1.0) / em;
        else if (m->usePwWeights == 2)  
            pwWeight[d]=v;

        MrBayesPrint("%s Pw Weight: %f \n", spacer, pwWeight[d]);

        /*  free allocations   */
        /*  helper matrices */

        ///MrBayesPrint("%s 1st chunk \n", spacer);        
        FreeSquareDoubleMatrix(V);
        FreeSquareDoubleMatrix(Vinv);
        FreeSquareComplexMatrix(Vc);
        FreeSquareComplexMatrix(Vcinv);
        free(tempPartitionPair);

        free(dw);
        free(iw);

        /*  counts */
        //MrBayesPrint("%s 2nd chunk \n", spacer);        
        free(counts);
        //MrBayesPrint("%s 2a \n", spacer);        
        for (i=0; i<(nSplits+1); i++)
            free(countIndex[i]);
        //MrBayesPrint("%s 2b \n", spacer);        
        free(countIndex);
        //MrBayesPrint("%s 2c \n", spacer);        
        //MrBayesPrint("%s 2d \n", spacer);        
        //free(n10);
        //MrBayesPrint("%s 2e \n", spacer);        
        //free(n11);
        //
        for (i=0; i<nPairs; i++) 
            free(tp1[i]);

        for (i=0; i<nPairs; i++) 
            free(tp[i]);

        free(tptemp);

        free(pwDists);
                      
        /*  hessian and jacobian */
        for (i=0; i<numBranches; i++) 
            { 
            free(H[i]);
            free(J[i]);
            free(HiJ[i]);
            free(Hinv[i]);
            }
        free(H);
        free(J);
        free(Hinv);
        free(HiJ);

        for (i=0; i<nSplits; i++) 
            { 
            free(D1L[i]);
            free(D1LP[i]);
            }
        free(D1L);
        free(D1LP);
        free(eigvals);
        free(eigvalsc);

        //MrBayesPrint("%s 6th free \n", spacer);
        for (i=0; i<numBranches; i++)
            free(PairBranch[i]);
        free(PairBranch);


        } /* end loop over numCurrentDivisions */
    return(0);
}
