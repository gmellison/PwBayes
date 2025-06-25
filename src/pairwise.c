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
#include "main.h"

#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#define LIKEPW_EPSILON        1.0e-300
#define LIKEPW_EPSILON        1.0e-300


// global variables for pairwise likelihoods:
// global variables for pw likelihoods:
MrBFlt nucFreqs[4];
MrBFlt doubletFreqs[4][4];

MrBFlt relRates[6];
int    tempNumStates;
int    numTriples;
int I=4;
int J=4;

int defFreqs = NO;
int defDoublets = NO;

MrBFlt alpha=0.1;
MrBFlt dist=0.1;

int *pairwiseCounts;
int **tripleCounts;
int numTrips;
int numPairs;

/* globals declared here */
int defPairwise=NO;
int defTriples=NO;
int allocPairwise=NO;
int allocTriples=NO;

/*  
int usePairwise=NO;
int useTriples=NO;
int useFullForAlpha=NO;
 */


// local prototypes:
int              toIdx(int x);
void             SetupQMat(MrBFlt **Q);
MrBFlt*          AllocateDists(int nPairs);
int***           AllocateDoubletCounts(int nPairs);
void             FreeDoubletCounts(int ***doubletCounts, int np);
int              CalcPairwiseDistsPolyTree(PolyTree *t, PairwiseDists *pd);
int              pairIdx(int i, int j, int n);
int              tripletIdx(int i, int j, int k, int n);
MrBFlt           CalcTripletLogLike_JC(PairwiseDists *pd);
int              FreeDists(MrBFlt* dists);
int              FreePairwiseDists(PairwiseDists* pd);
double           CalcPairwiseLogLike_Full(PairwiseDists *pd);
int              FreeTriples(void);


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

int tripletIdx(int i, int j, int k, int n) {

    int l;

    if (i == j || j == k || i == k || i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n || n < 3) {
        MrBayesPrint("Error in triplet indexing; i=%d, j=%d,k=%d,n=%d \n", i,j,k,n);
        return(-1);
    }

    /*  if i > j, swap i,j */
    if (i > j) {
        l = i;
        i = j;
        j = l;
    } 

    /*  if i > k, swap i,k so i is for sure the minimum */
    if (i > k) {
        l = i;
        i = k;
        k = l;
    } 

    /*  now if j > k, swap j,k */
    if (j > k) {
        l = j;
        j = k;
        k = l;
    } 

    return(n*(n-1)*(n-2)/6 - (n-i-1)*(n-i-2)*(n-i-3)/6 - (n-j-1)*(n-j-2)/2 - (n-k-1) - 1) ;
}

int triplePos(int i, int j, int k) {
    return i*16 +j*4 + k;
}


/*  *****************
 *
 *
 *  Pw Rate Estimation & standalone likelihood calculations
 *
 *
 *  *****************  */

int DoEstQPairwise(void) {

    MrBayesPrint("%s Running 'Estqpairwise' command \n", spacer);

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }
           
    if (defFreqs == NO) 
    {
        MrBayesPrint("%s Getting Nucleotide Frequences \n", spacer);
        CountFreqs(); // sets nucleotide frequencies in global object 'nucFreqs'
    }

    if (defDoublets == NO) 
    {
        MrBayesPrint("%s Getting Doublet Counts \n", spacer);
        CountDoublets(numTaxa);
    }

    MrBayesPrint("%s Estimating Q Pairwise  \n", spacer);

    int i,j,k;
    int t1, t2;
    int r1, r2, ridx;
    MrBFlt tau;
    MrBFlt *ratesEst;

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

    int numPairs = numTaxa*(numTaxa-1)/2;
    ratesEst = malloc(6*numPairs*sizeof(MrBFlt));

    for (t1=0; t1<(numTaxa-1); t1++) {
        for (t2=t1+1; t2<numTaxa; t2++) {
 
            // use mb machinery to compute eigens
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

            k = pairIdx(t1,t2,numTaxa);

            // set up matrix of empirical transitions 
            //   probabilities, 'Q' (will be modified in place)  
            for (i=0;i<4;i++) {
                for (j=0;j<4;j++){
                    F[i][j] = 1.0 * pairwiseCounts[tIdx(k,i,j,I,J)] / (2*numChar*nucFreqs[j]);
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


            int isComplex=GetEigens(4,F,la,laC,V,Vinv,Vc,Vcinv);
            (void)isComplex;
            
            // diagonal matrix of log^{lambda_i}, lambdas are eigvals of Q
            for (i=0;i<4;i++) {
                LaLog[i][i] = log(la[i]);
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
                Q[i][i]=-1*rowsum;
            }
        
            for (i=0;i<4;i++)
                tau += -1.0 * Q[i][i] * nucFreqs[i];
        
            MultiplyMatrixByScalar(4, Q, 1.0/tau, Q);
        
            MrBayesPrint("%s Tau Hat (%d,%d): %f \n", spacer, t1,t2, tau);
        
            // get the individual rates: 
            for (r1=1; r1<4; r1++) {
                for (r2=0; r2<r1; r2++) {
                    ridx = pairIdx(r1,r2,4);
                    ratesEst[dIdx(k,ridx,6)] = Q[r1][r2];
                }
            }

            MrBayesPrint("%s Estimated Relative Rates (pair %d): ", spacer, k);
            for (i=0;i<6;i++)
                MrBayesPrint(" %f", ratesEst[dIdx(k,i,6)]);
            MrBayesPrint("\n");

            // free all matrices 
            FreeSquareDoubleMatrix(V);
            FreeSquareDoubleMatrix(Vinv);
            FreeSquareDoubleMatrix(Q);
            FreeSquareComplexMatrix(Vc);
            FreeSquareComplexMatrix(Vcinv);
            FreeSquareDoubleMatrix(Temp);
            FreeSquareDoubleMatrix(F);
        }
    }


    /*  average the rates across pairs */
    MrBFlt r6;
    MrBFlt ratesOut[6] = {0.0}; 

    
    /*  Normalize the rates:  */
    for (k=0; k<numPairs; k++) {

        r6=ratesEst[dIdx(k,5,6)];
        for (i=0;i<6;i++) {
            ratesEst[dIdx(k,i,6)]/=r6;
            ratesOut[i]+=ratesEst[dIdx(k,i,6)]/(1.0*numPairs);
        }
    }

    MrBayesPrint("%s Estimated Relative Rates (Averaged): ", spacer);
    for (i=0;i<6;i++)
        MrBayesPrint(" %f", ratesOut[i]);
    MrBayesPrint("\n");

    free(ratesEst);

    return(NO_ERROR);
}



double CalcPairwiseLogLike_Reduced(PairwiseDists *pd) {

    int i,j,k;
    int ri;
    MrBFlt tp;
    MrBFlt ll=0.0;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    /* Compute average pairwise dist:  */
    MrBFlt dist = 0.0;
    for (k=0;k<pd->nPairs;k++) 
        dist += (pd->dists)[k];

    dist = dist/(1.0*pd->nPairs);


    /*  pool doublet counts:  */
    int **doubletCountsPooled;
    doubletCountsPooled = AllocateSquareIntegerMatrix(4);

    for (k=0; k<pd->nPairs; k++) {
        for (i=0; i<4; i++) 
            for (j=0; j<4; j++) 
                doubletCountsPooled[i][j] += pairwiseCounts[tIdx(k,i,j,I,J)];
    }
  
    // compute gtr transition probabilities
    // set up matrices for taking the matrix exponent
    MrBFlt **V         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Vinv      = AllocateSquareDoubleMatrix(4);
    MrBFlt **LaExp     = AllocateSquareDoubleMatrix(4);
    MrBFlt **TempMat   = AllocateSquareDoubleMatrix(4);

    MrBComplex **Vc    = AllocateSquareComplexMatrix(4);    // eigenvectors
    MrBComplex **Vcinv = AllocateSquareComplexMatrix(4);    // inverse eigvect matrix
    MrBFlt *la         = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // eigenvalues
    MrBFlt *laC        = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // 

    MrBFlt **Q         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Qtausr    = AllocateSquareDoubleMatrix(4);

    // allocate array of 4x4 matrices -- one for each rate category
    MrBFlt **pTrans[100];
    for (int ri=0; ri<4; ri++)
            pTrans[ri]=AllocateSquareDoubleMatrix(4);

    SetupQMat(Q);

    // calculate the site pattern probability for each 
    //   dg4 rate category
    for (ri=0; ri<4; ri++) {
    
        // probability transition matrix for site rate i:
        MultiplyMatrixByScalar(4, Q, dist*dg4[ri], Qtausr);  

        int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
        if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!! \n");
           
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
        MultiplyMatrices(4,TempMat,Vinv,pTrans[ri]);
    }

 
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++) 
            MrBayesPrint("%s Q(%d,%d) %f \n", spacer, i,j, pTrans[1][i][j]);
    
    // now sum over all the doublet patterns
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
         
            // reset tp, and sum up the conditional 
            //   trans probs over discrete rate categories
            tp=0.0; 
            for (int ri=0; ri<4; ri++) {
                tp += 0.25 * pTrans[ri][i][j];
            }

            // calculate triple probability by summing over 
            //   conditional probs given rate categories 
            ll += doubletCountsPooled[i][j] * log(tp);
            
        }
    }

    FreeSquareDoubleMatrix(TempMat);
    for (int ri=0;ri<4;ri++)
        FreeSquareDoubleMatrix(pTrans[ri]);


    return ll;

}

double CalcPairwiseLogLike_Full(PairwiseDists *pd) {

    int i,j;
    int ri;
    int t1, t2, k;
    MrBFlt tp;
    MrBFlt ll=0.0;
    MrBFlt dist;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    MrBFlt **Q = AllocateSquareDoubleMatrix(4);
    SetupQMat(Q);

    for (t1=1; t1<pd->nTaxa; t1++)
        {
        for (t2=0; t2<t1; t2++)
            {

            k = pairIdx(t1,t2,pd->nTaxa);    
            dist = pd->dists[k];

            // compute gtr transition probabilities
            // set up matrices for taking the matrix exponent
            MrBFlt **V         = AllocateSquareDoubleMatrix(4);
            MrBFlt **Vinv      = AllocateSquareDoubleMatrix(4);
            MrBFlt **LaExp     = AllocateSquareDoubleMatrix(4);
            MrBFlt **TempMat   = AllocateSquareDoubleMatrix(4);

            MrBComplex **Vc    = AllocateSquareComplexMatrix(4);    // eigenvectors
            MrBComplex **Vcinv = AllocateSquareComplexMatrix(4);    // inverse eigvect matrix
            MrBFlt *la         = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // eigenvalues
            MrBFlt *laC        = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // 

            MrBFlt **Qtausr    = AllocateSquareDoubleMatrix(4);

            // allocate array of 4x4 matrices -- one for each rate category
            MrBFlt **pTrans[100];
            for (int ri=0; ri<4; ri++)
                    pTrans[ri]=AllocateSquareDoubleMatrix(4);

            // calculate the site pattern probability for each 
            //   dg4 rate category
            for (ri=0; ri<4; ri++) {
            
                // probability transition matrix for site rate i:
                MultiplyMatrixByScalar(4, Q, dist*dg4[ri], Qtausr);  

                int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
                if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!! \n");
                   
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
                MultiplyMatrices(4,TempMat,Vinv,pTrans[ri]);
            }
       
            // now sum over all the doublet patterns
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                 
                    // reset tp, and sum up the conditional 
                    //   trans probs over discrete rate categories
                    tp=0.0; 
                    for (int ri=0; ri<4; ri++) {
                        tp += 0.25 * pTrans[ri][i][j];
                    }

                    // calculate triple probability by summing over 
                    //   conditional probs given rate categories 
                    MrBayesPrint("%s k=%d, i=%d, j=%d, tIdx: %d, Count: %d, log(p): %f \n", 
                            spacer,
                            k,i,j,
                            tIdx(k,i,j,I,J), 
                            pairwiseCounts[tIdx(k,i,j,I,J)],  
                            log(tp));
                    ll += pairwiseCounts[tIdx(k,i,j,I,J)] * log(tp);
                    
                }

            }

            FreeSquareDoubleMatrix(TempMat);
            for (int ri=0;ri<4;ri++)
                FreeSquareDoubleMatrix(pTrans[ri]);
            }
        }


    return ll;

}

int DoPairwiseLogLike(void) {

    MrBayesPrint("Running 'Pairwiseloglike' command.\n");

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }

    if (defFreqs==NO) 
    {
        CountFreqs();
    }

    if (defDoublets==NO) 
    {
        CountDoublets(numTaxa);
        MrBayesPrint("Done Counting Doublets \n");
    }

    if (numUserTrees==0) {
        MrBayesPrint ("%s   No user trees defined! \n", spacer);
        return (ERROR);
    }

    int l;
    MrBFlt ll;
    PairwiseDists *pd; 

    MrBayesPrint ("%s   Calculating Pw likelihood for %d defined user trees. \n", spacer, numUserTrees);

    for (l=0; l<numUserTrees; l++) {
        pd = AllocatePairwiseDists();
        InitPairwiseDistsPolyTree(userTree[l],pd);

        ll=CalcPairwiseLogLike_Full(pd);
        MrBayesPrint("Pairwise Log-Likelihood: %f \n", ll);
        FreePairwiseDists(pd);
    }


    return(NO_ERROR);

}

/*  *****************
 *
 *
 *  Triplet composite likelihood functions 
 *
 *
 *  *****************  */

                                  
int FreeTriples(void) 
{
    int i;

    for (i=0; i<numTrips; i++)
        if (tripleCounts[i] != NULL)
            free(tripleCounts[i]);
        else 
            return (ERROR);

    if (tripleCounts != NULL)
        free(tripleCounts);
    else 
        return (ERROR);
    
    return (NO_ERROR);
}




int CountTriplets(void) {

    int             i,j,k,c,id1,id2,id3,tripIdx;
    int             *counts;
    
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

    if (defMatrix == NO) 
        {
        MrBayesPrint("%s Matrix needs to be defined before counting triplets  \n", spacer);
        return(ERROR);
        }  
    
    if (numDefinedPartitions > 1)    
        { 
        MrBayesPrint("%s Pairwise count likelihood only implemented for a single partition  \n", spacer);
        return(ERROR);
        }

    if (memAllocs[ALLOC_TRIPLES] == YES) 
        {
        FreeTriples();
        tripleCounts=NULL;
        memAllocs[ALLOC_TRIPLES] = NO;
        }


    numTrips=(int)numTaxa*(numTaxa-1)*(numTaxa-2)/6;
    MrBayesPrint("NumTriplets: %d \n", numTrips);

    /* first allocate pairwise doublet counts:  */
    tripleCounts=(int**)SafeMalloc(numTrips * sizeof(int*));
    for (i=0; i<numTrips; i++)
            tripleCounts[i]=(int*)SafeMalloc(64*sizeof(int));

    if (!pairwiseCounts)
        {
        MrBayesPrint("%s Problem allocating pairwise counts! \n", spacer);
        free(pairwiseCounts);
        return(ERROR);
        }

    /* now count the doublets across taxa pairs */
    tripIdx=0;
    for (i=0; i<(numTaxa-2); i++)
        {
        for (j=i+1; j<(numTaxa-1); j++)
            {
            for (k=j+1; k<numTaxa; k++)
                {
                counts=tripleCounts[tripIdx++];
                for (c=0;c<numChar;c++)
                    {
                    if (matrix[pos(i,c,numChar)]==GAP | matrix[pos(j,c,numChar)]==GAP | matrix[pos(k,c,numChar)]==GAP)
                        continue;

                    /* nucleotides at position x of sequences i & j   */        
                    id1=toIdx(matrix[pos(i,c,numChar)]);
                    id2=toIdx(matrix[pos(j,c,numChar)]);
                    id3=toIdx(matrix[pos(k,c,numChar)]);

                    /* increment the count of that nucleotide pair, at pair k  */
                    counts[tIdx(id1,id2,id3,4,4)]++;
                    }
                }
            }
        }

    return (NO_ERROR);
}


int CalcTripletCnDists(int division, int chain)
{
    int         i,j,k,tripIdx;
    MrBFlt      *pwDists, *cnDists, pwd1, pwd2, pwd3; 
    ModelInfo   *m;  

    m=&modelSettings[division];
    pwDists=m->pwDists[chain];
    cnDists=m->tripleCnDists[chain];

    tripIdx=0;
    for (i=0; i<(numTaxa-2); i++)
        {
        for (j=i+1; j<(numTaxa-1); j++)
            {
            pwd1=pwDists[pairIdx(i,j,numTaxa)];
            for (k=j+1; k<numTaxa; k++)
                {
                pwd2=pwDists[pairIdx(i,k,numTaxa)];
                pwd3=pwDists[pairIdx(j,k,numTaxa)];

                cnDists[tripIdx++]=(pwd1+pwd2-pwd3)/(2.0);
                cnDists[tripIdx++]=(pwd1+pwd3-pwd2)/(2.0);
                cnDists[tripIdx++]=(pwd2+pwd3-pwd1)/(2.0);
                }
            }
        }

    
    return(NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   TiProbsTriplet_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsTriplet_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    int         i, j, k, l, p, index, tripIdx;
    MrBFlt      *tripleDists; 
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange;
    CLFlt       *tiP;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    /* find pw dists and transition probabilities */
    tripleDists = m->tripleCnDists[chain];

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

    tripIdx=0;
    for (p=0; p<numTrips; p++)
        {

        for (l=0;l<3;l++)
            {
            
            length = tripleDists[tripIdx];
            tiP = m->tripleTiProbs[m->tripDistIndex[chain][tripIdx]];
            tripIdx++;

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
        }

    return (NO_ERROR);
}

/*-----------------------------------------------------------------
|
|   TiProbsTriplet_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TripletProbs_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    int         k, p, w, spIdx, tripLengthIdx, x,y,z, idx1, idx2, idx3;
    CLFlt       numCats;
    CLFlt       *tiP1, *tiP2, *tiP3, *tripProb;   
    ModelInfo   *m;
    int         indexStep, indexStart;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    numCats=((CLFlt)m->numRateCats);

    /*  tripProbTemp=(CLFlt*)SafeMalloc(4 * 4 * 4 * sizeof(CLFlt)); */

    indexStep=4*4;

    tripLengthIdx=0;
    for (p=0; p<numTrips; p++)
        {
        indexStart=0;
        tripProb=m->tripleProbs[m->tripIndex[chain][p]];

        /*  initialize triplet probabilities to 0: */
        for (k=0; k<m->tripleProbsLength; k++)
                tripProb[k]=0.0;

        tiP1 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
        tiP2 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
        tiP3 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
 
        for (k=0; k<m->numRateCats; k++)
            {

            spIdx=0;
            idx1=indexStart;   
            for (x=0; x<4; x++)
                {
                idx2=indexStart; 
                for (y=0; y<4; y++)
                    {
                    idx3=indexStart;
                    for (z=0; z<4; z++)
                        {
                        for (w=0; w<4; w++)
                            {
                                /*  internal triplet nucleotide is w */
                                /*  need transition probabilities w->x,w->y,w->z  */
                                tripProb[spIdx] += 0.25 * tiP1[idx1+w] * tiP2[idx2+w] * tiP3[idx3+w] / numCats;
                                /*
                                if (chain ==0 && x == 0 && y == 0 && z == 2 && p == 0) {

                                        MrBayesPrint("Site Pattern AAC|w=%d, triplet 1: ---- \n ", w);
                                        MrBayesPrint("  Triple Probs | rate cat = %d:\n ", k);
                                        MrBayesPrint("  tiP1[idx1+w]: %f,  tiP2[idx2+w]: %f tiP3[idr3+w]: %f  \n ", tiP1[idx1+w], tiP2[idx2+w], tiP3[idx3+w] );

                                        dummy += 0.25 * tiP1[idx1+w] * tiP2[idx2+w] * tiP3[idx3+w] / numCats ;
                                        if (w == 3) MrBayesPrint("  triplet prob for rate cat %d: %f \n\n", k, dummy);
                                }  */
                            }
                        spIdx++;
                        idx3+=4;
                        } /*  end of z loop */
                    idx2+=4;
                    } /*  end of y loop  */
                idx1+=4;
                } /*  end of x loop */

            indexStart+=indexStep;
            } /* end rate category loop  */

        /*   tripProbTemp[tripIdx]=;   */
        } /*  end loop over triplets */

    /*  free(tripProbTemp);  */
    return (NO_ERROR);
}

int Likelihood_Triples (int division, int chain, MrBFlt *lnL)
{
    /* calculate Jukes Cantor likelihood, using pw counts and tps. */
    
    int         i, j, k, idx, p, nijk;
    MrBFlt      like;
    CLFlt       *tripP, pijk;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    (*lnL) = 0.0;

    for (p=0; p<numTrips; p++)
        {

        /* find transition probabilities */
        tripP = m->tripleProbs[m->tripIndex[chain][p]];

        idx=0;
        for (i=0; i<4; i++) 
            {
            for (j=0; j<4; j++) 
                {
                for (k=0; k<4; k++)
                    {
                    like = 0.0;
                    nijk=tripleCounts[p][idx];
                    pijk=tripP[idx++];

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

                    *lnL+=like;
                    }
                }
            }
        }

    return (NO_ERROR);
}



MrBFlt LogLikeTriplet_Alpha(int chain)
{
    MrBFlt     lnL;
    int d=0; /*  for now, only implemented for a single DNA partition */

    CalcTripletCnDists(d, chain); 
    TiProbsTriplet_JukesCantor(d, chain); 
    TripletProbs_JukesCantor(d,chain);
    TIME(Likelihood_Triples(d,chain,&lnL),CPULilklihood); 
    return(lnL);
}

int DoTripletLogLike(void) {

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }

    if (numUserTrees == 0)
    {
        MrBayesPrint ("%s   No user trees defined. \n", spacer);
        return (ERROR);
    }

    if (defFreqs==NO) 
    {
        CountFreqs();
    }

    if (defTriples==NO)
    {
        CountTriplets();
    }

    PairwiseDists* pd;
    int ti;
    MrBFlt ll=0.0;

    for (ti=0; ti<numUserTrees; ti++) {
        pd = AllocatePairwiseDists();
        InitPairwiseDistsPolyTree(userTree[ti],pd);

        ll = CalcTripletLogLike_JC(pd);
        MrBayesPrint("Triplet Quasi Log-Likelihood for Tree %d: %f", ti, ll);
        FreePairwiseDists(pd);
    }

    return(NO_ERROR);

}


MrBFlt CalcTripletLogLike_JC(PairwiseDists *pd) {

    double al=0.5;
    int i,j,k,ri;
    int seq1, seq2, seq3;
    int pairidx, tripidx, nucidx;
     

    MrBFlt dg4[4];
    DiscreteGamma(dg4,al,al,4,1);

    MrBFlt pii[3][4];
    MrBFlt pij[3][4];

    MrBFlt tripleProbs[64];
    MrBFlt transitionProbs[3][64];

    MrBFlt pwdist[3];
    MrBFlt cndist[3];
 
    MrBFlt probByRateCat[4];
    MrBFlt ll=0.0;
 
    /* 
     * loop over alignment sequence triples 
     */ 
    for (seq1=0; seq1<(pd->nTaxa-2); seq1++) {
        for (seq2=seq1+1; seq2<(pd->nTaxa-1); seq2++) {
            for (seq3=seq2+1; seq3<pd->nTaxa; seq3++) {

                tripidx = tripletIdx(seq1,seq2,seq3,pd->nTaxa);

                /*
                 * get the distances between triple taxa (pw and cn)
                 */
                pwdist[0]=pd->dists[pairIdx(seq1,seq2,pd->nTaxa)];
                pwdist[1]=pd->dists[pairIdx(seq1,seq3,pd->nTaxa)];
                pwdist[2]=pd->dists[pairIdx(seq2,seq3,pd->nTaxa)];
                MrBayesPrint("Triple Dists: %f, %f, %f \n", pwdist[0], pwdist[1], pwdist[2]);

                cndist[0] = (pwdist[0] + pwdist[1] - pwdist[2])/2.0;
                cndist[1] = (pwdist[0] + pwdist[2] - pwdist[1])/2.0;
                cndist[2] = (pwdist[1] + pwdist[2] - pwdist[0])/2.0;
                MrBayesPrint("Centernode Dists: %f, %f, %f \n", cndist[0], cndist[1], cndist[2]);

                /*  
                 *  calc JC transition probs for each centernode dist
                 */
                for (ri=0; ri<4; ri++) {
                    for (pairidx=0; pairidx<3; pairidx++) {
                        MrBFlt etr = exp(-(4.0/3) * cndist[pairidx] * dg4[ri]);
                        pii[pairidx][ri]=(1.0/4) + (3.0/4) * etr;
                        pij[pairidx][ri]=(1.0/4) - (1.0/4) * etr;

                        /*
                        MrBayesPrint("  Pij[%d][%d]: %f", pairidx, ri, pij[pairidx][ri]);
                        MrBayesPrint("  Pii[%d][%d]: %f", pairidx, ri, pii[pairidx][ri]);
                         */

                        // set up array of nucleotide transition probabilities                                
                        for (i=0;i<4;i++) {
                            for (j=0;j<4;j++) {
                                if (i==j) {
                                    transitionProbs[pairidx][tIdx(i,j,ri,4,4)]=pii[pairidx][ri];
                                } else {
                                    transitionProbs[pairidx][tIdx(i,j,ri,4,4)]=pij[pairidx][ri];
                                }
                            }
                        }
                    }
                }

                /*  
                *  calc JC transition probs for each centernode dist
                */
                for (i=0; i<4; i++) {
                    for (j=0; j<4; j++) {
                        for (k=0; k<4; k++) {

                            nucidx = tIdx(i,j,k,4,4);

                            // no need to calculate probs for not present in data
                            if (tripleCounts[tripidx][nucidx] == 0) continue;  

                            // ensure triple site pattern prob is init to 0:
                            tripleProbs[nucidx] = 0.0;

                            // calculate the site pattern probability for each 
                            //   dg4 rate category
                            for (ri=0; ri<4; ri++) {

                                // reset helper rate cat array at current index
                                probByRateCat[ri]=0.0;

                                // sum over the 4 possible internal nucleotides, 
                                // given rate cat ri:   
                                for (int nucx=0; nucx<4; nucx++) {
                                        probByRateCat[ri]+=
                                                    nucFreqs[nucx]
                                                    *transitionProbs[0][tIdx(nucx,i,ri,4,4)]
                                                    *transitionProbs[1][tIdx(nucx,j,ri,4,4)]
                                                    *transitionProbs[2][tIdx(nucx,k,ri,4,4)];
                                }

                                // calculate triple probability by summing over 
                                //   conditional probs given rate categories 
                                tripleProbs[nucidx]+=0.25*probByRateCat[ri];
                            }

                            ll += tripleCounts[tripidx][nucidx] * log(tripleProbs[nucidx]);
                        }
                    }
                }
            }
        }
    }
    return(ll);
}

  
MrBFlt CalcTripletLogLike_JC_Reduced(PairwiseDists *pd, MrBFlt alpha) {


    MrBFlt tripleProbs[64], transitionProbs[64];

    int i,j,k, ri;
    int ti;
    int spIdx;

    MrBFlt dist;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    MrBFlt pii[4];
    MrBFlt pij[4];

    /*  Calculate avg dist: 
     */
    dist = 0.0;
    for (k=0;k<pd->nPairs;k++) 
        dist += (pd->dists)[k];

    dist = dist/(2.0*pd->nPairs);


    /*  pool triplet counts:  */
    int tripleCountsPooled[64];

    for (ti=0; ti<pd->nPairs; ti++)
        for (i=0; i<4; i++) 
            for (j=0; j<4; j++) 
                for (k=0;k<4;k++) 
                    tripleCountsPooled[tIdx(i,j,k,4,4)] += tripleCounts[ti][tIdx(i,j,k,4,4)];
    
 
    for (int i=0; i<4; i++) {
        MrBFlt etr = exp(-(4.0/3) * dist* dg4[i]);
        pii[i]=(1.0/4) + (3.0/4) * etr;
        pij[i]=(1.0/4) - (1.0/4) * etr;
    }

    // set up array of transition probabilities
    for (i=0;i<4;i++) {
            for (j=0;j<4;j++) {
                    for (int ri=0;ri<4;ri++) {
                            if (i==j) {
                                transitionProbs[triplePos(i,j,ri)]=pii[ri];
                            } else {
                                transitionProbs[triplePos(i,j,ri)]=pij[ri];
                            }
                    }
            }
    }

    MrBFlt probsByRateCat[4];
    MrBFlt ll=0.0;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {

                spIdx=tIdx(i,j,k,4,4);

                // no need to calculate probs for not present in data
                if (tripleCounts[spIdx] == 0) continue;  

                // ensure triple site pattern prob is init to 0:
                tripleProbs[spIdx] = 0.0;

                // calculate the site pattern probability for each 
                //   dg4 rate category
                for (ri=0; ri<4; ri++) {

                    // reset helper rate cat array at current index
                    probsByRateCat[ri]=0.0;

                    // sum over the 4 possible internal nucleotides, 
                    // given rate cat ri:   
                    for (int nucx=0; nucx<4; nucx++) {
                            probsByRateCat[ri]+=
                                        nucFreqs[nucx]
                                        *transitionProbs[triplePos(nucx,i,ri)]
                                        *transitionProbs[triplePos(nucx,j,ri)]
                                        *transitionProbs[triplePos(nucx,k,ri)];
                    }

                    // calculate triple probability by summing over 
                    //   conditional probs given rate categories 
                    tripleProbs[spIdx]+=0.25*probsByRateCat[ri];
                }

                ll += tripleCountsPooled[spIdx] * log(tripleProbs[spIdx]);
            }
        }
    }
    return (ll);
}





/*  *****************
 *
 *
 *  Pw likelihood for Mrb MCMC
 *   
 *
 *  *****************  */

/* Inits a PairwiseDists struct from a polytree  
 *
 *  
 */

PairwiseDists* AllocatePairwiseDists(void) {
    PairwiseDists *pd;
    pd=(PairwiseDists*)SafeCalloc(1,sizeof(PairwiseDists));
    return(pd);
}

/*
void InitPairwiseDists(Tree *t, PairwiseDists *pd) 
{
    pd->nTaxa=(t->nNodes-t->nIntNodes);
    pd->nPairs=(pd->nTaxa)*(pd->nTaxa-1)/2;
    pd->dists=AllocateDists(pd->nTaxa);
    CalcPairwiseDists_ReverseDownpass(t, pd);
}
*/

void InitPairwiseDistsPolyTree(PolyTree *pt, PairwiseDists *pd) 
{
    pd->nTaxa=(pt->nNodes-pt->nIntNodes);
    pd->nPairs=(pd->nTaxa)*(pd->nTaxa-1)/2;
    pd->dists=AllocateDists(pd->nTaxa);
    CalcPairwiseDistsPolyTree(pt, pd);
}

int FreePairwiseDists(PairwiseDists* pd) 
{
    FreeDists(pd->dists);
    free(pd);
    return(NO_ERROR);
}

//int InitPairwiseCondLikes(int numLocalChains) 
//{
//    int         d, i, j, c, k, id1, id2,
//                indexStep, pwIdx, tripIdx;
//    ModelInfo   *m;
//
//    for (d=0; d<numCurrentDivisions; d++)
//        {
//        m = &modelSettings[d];
//        if (m->usePairwise == NO) 
//            continue;
//        
//        if (m->dataType != DNA)
//            {
//            MrBayesPrint("%s attempt to init pairwise for division with non DNA data", spacer);
//            return ERROR; 
//            }
//
//        /* find size of tree */
//        //nIntNodes = GetTree(m->brlens, 0, 0)->nIntNodes;
//        //nNodes = GetTree(m->brlens, 0, 0)->nNodes;
//
//        m->numPairs = (numLocalTaxa) * (numLocalTaxa - 1) / 2;
//        m->tiProbsPwLength = m->numModelStates * m->numModelStates * m->numTiCats;
//        m->doubletProbsLength = m->numModelStates * m->numModelStates;
//        //m->numTiProbs = (numLocalChains + 1) * nNodes;
//        m->numTiProbsPw = (numLocalChains + 1) * m->numPairs;
//        m->numDoubletProbs = (numLocalChains + 1) * m->numPairs;
//        
//        /*  allocate space for pw distances */
//        m->pwDists = (MrBFlt**) SafeMalloc(numLocalChains * sizeof(MrBFlt*));
//        for (i=0; i<numLocalChains; i++)
//            m->pwDists[i] = (MrBFlt*) SafeMalloc(m->numPairs * sizeof(MrBFlt));
//
//        /*  allocate space for pw ti probs  */
//        m->tiProbsPw = (CLFlt**) SafeMalloc(m->numTiProbsPw * sizeof(CLFlt*));
//        if (!m->tiProbs)
//            return (ERROR);
//        for (i=0; i<m->numTiProbsPw; i++)
//            {
//            m->tiProbsPw[i] = (CLFlt*)SafeMalloc(m->tiProbsPwLength * sizeof(CLFlt));
//            if (m->tiProbsPw[i] == NULL)
//                 return (ERROR);
//            }
//
//        /*  allocate space for pw doublet probs  */
//        m->doubletProbs = (CLFlt**) SafeMalloc(m->numDoubletProbs * sizeof(CLFlt*));
//        if (!m->doubletProbs)
//            return (ERROR);
//        for (i=0; i<m->numDoubletProbs; i++)
//            {
//            m->doubletProbs[i] = (CLFlt*)SafeMalloc(m->doubletProbsLength * sizeof(CLFlt));
//            if (!m->doubletProbs[i])
//                return (ERROR);
//            }
//
//        /* allocate and set indices from chain/pair to pw probs */
//        m->pwIndex = (int **) SafeMalloc (numLocalChains * sizeof(int *));
//        if (!m->pwIndex)
//            return (ERROR);
//        for (i=0; i<numLocalChains; i++)
//            {
//            m->pwIndex[i] = (int *) SafeMalloc (m->numPairs * sizeof(int));
//            if (!m->pwIndex[i])
//                return (ERROR);
//            }
//
//        indexStep=1;
//        /* set up pw indices */
//        pwIdx = 0;
//        for (i=0; i<numLocalChains; i++)
//            {
//            for (j=0; j<m->numPairs; j++)
//                {
//                m->pwIndex[i][j] = pwIdx;
//                pwIdx += indexStep;
//                }
//            }
//
//        m->pwCounts=(int*)SafeMalloc(m->numPairs * 16 * sizeof(int));
//        if (!m->pwCounts)
//            {
//            MrBayesPrint("%s Problem allocating pairwise counts! \n", spacer);
//            free(m->pwCounts);
//            return(ERROR);
//            }
//
//        ///memAllocs[ALLOC_PAIRWISE]=YES;
//
//        /* now count the doublets across taxa pairs */
//        for (i=0; i<(numTaxa-1); i++)
//            {
//            for (j=i+1; j<numTaxa; j++)
//                {
//                k=pairIdx(i,j,numTaxa);
//                for (c=0;c<numChar;c++)
//                    {
//                    if (matrix[pos(i,c,numChar)]==GAP | matrix[pos(j,c,numChar)]==GAP)
//                        continue;
//
//                    if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum]!=d+1) 
//                        continue;
//
//                    /* nucleotides at position x of sequences i & j   */        
//                    id2=toIdx(matrix[pos(j,c,numChar)]);
//                    id1=toIdx(matrix[pos(i,c,numChar)]);
//
//                    /* increment the count of that nucleotide pair, at pair k  */
//                    m->pwCounts[tIdx(k,id1,id2,4,4)]++;
//                    }
//                }
//            }
//
//        memAllocs[ALLOC_PAIRWISE] = YES;
//
//        /*  allocate triplet stuff if necessary */
//        if (!m->useTriples) 
//            continue;
//
//        m->numTrips = (numLocalTaxa) * (numLocalTaxa - 1) * (numLocalTaxa - 2)/ 6;
//        m->tiProbsTripLength = m->numModelStates * m->numModelStates * m->numTiCats;
//        m->tripleProbsLength = m->numModelStates * m->numModelStates * m->numModelStates;
//        m->numTiProbsTrip = (numLocalChains + 1) * m->numTrips * 3;
//        m->numTripleProbs = (numLocalChains + 1) * m->numTrips;
//    
//        /*  allocate triple cn distances */
//        m->tripleCnDists=(MrBFlt**)SafeMalloc(numLocalChains * sizeof(MrBFlt*));
//        for (i=0;i<numLocalChains;i++)
//            m->tripleCnDists[i]=(MrBFlt*)SafeMalloc(m->numTrips * 3 * sizeof(MrBFlt));
//
//        /*  allocate triple ti probabilities */
//        m->tripleTiProbs=(CLFlt**)SafeMalloc(m->numTiProbsTrip * sizeof(CLFlt*));
//        for (i=0;i<m->numTiProbsTrip;i++)
//            m->tripleTiProbs[i]=(CLFlt*)SafeMalloc(m->tiProbsTripLength * sizeof(CLFlt));
//        MrBayesPrint("TripleTiProbs Size: %d, %d \n", m->numTiProbsTrip, m->tiProbsTripLength);
//
//        /*  allocate triple site pattern probabilities */
//        m->tripleProbs=(CLFlt**)SafeMalloc(m->numTripleProbs * sizeof(CLFlt*));
//        for (i=0;i<m->numTripleProbs;i++)
//            m->tripleProbs[i]=(CLFlt*)SafeMalloc(m->tripleProbsLength * sizeof(CLFlt));
//
//         /*  allocate triple indices  */
//        m->tripIndex = (int **) SafeMalloc (numLocalChains * sizeof(int *));
//        if (!m->tripIndex)
//            return (ERROR);
//        for (i=0; i<numLocalChains; i++)
//            {
//            m->tripIndex[i] = (int *) SafeMalloc (m->numTrips * sizeof(int));
//            if (!m->tripIndex[i])
//                return (ERROR);
//            }
//
//        /* set up triple indices */
//        tripIdx = 0;
//        for (i=0; i<numLocalChains; i++)
//            {
//            for (j=0; j<m->numTrips; j++)
//                {
//                m->tripIndex[i][j] = tripIdx;
//                tripIdx += indexStep;
//                }
//            }
//            
//      
//        /* set pwWeight temporarily. It'll be updated once the chain has run for some time  */ 
//        if (m->usePwWeights)  
//            m->pwWeight=1.0;  
//
//        }
//    return NO_ERROR;
//}

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

        /*  free triplet dists and probabilities  */
        if (m->useTriples) 
            {
            if (m->tripleCnDists)
                {
                for (j=0; j<numLocalChains; j++)
                    {
                    if (m->tripleCnDists[j])
                        free(m->tripleCnDists[j]);
                    }
                free(m->tripleCnDists);
                }

            if (m->tripleTiProbs)
                {
                for (j=0; j<numLocalChains; j++)
                    if (m->tripleTiProbs[j]) 
                        free(m->tripleTiProbs[j]);
                free(m->tripleTiProbs);
                }

            if (m->tripleProbs)
                {
                for (j=0; j<numLocalChains; j++)
                    if (m->tripleProbs[j]) 
                        free(m->tripleProbs[j]);
                free(m->tripleProbs);
                }

            if (m->tripIndex)
                {
                for (j=0; j<numLocalChains; j++)
                    {
                    if (m->tripIndex[j]) /*  chain pw probs */
                        free(m->tripIndex[j]);
                    }
                free(m->tripIndex);
                }

            if (m->tripDistIndex)
                {
                for (j=0; j<numLocalChains; j++)
                    {
                    if (m->tripDistIndex[j]) /*  chain pw probs */
                        free(m->tripDistIndex[j]);
                    }
                free(m->tripDistIndex);
                }
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



int CalcPairwiseDistsPolyTree(PolyTree *t, PairwiseDists *pd)
{
    int         i,j,a,d,k;
    PolyNode    *p;
    double      x;

    int numExtNodes = t->nNodes - t->nIntNodes;

    /*  We'll calculate (in distsTemp) all node dists including internal nodes  */
    MrBFlt **distsTemp;
    distsTemp=(MrBFlt**)malloc(t->nNodes*sizeof(MrBFlt*)); 
    for (k=0;k<t->nNodes;k++)
            distsTemp[k]=(MrBFlt*)malloc(t->nNodes*sizeof(MrBFlt));

    /*  make sure dists are init to 0  */
    for (i=0; i<t->nNodes; i++) 
        for (j=0; j<t->nNodes; j++)
                distsTemp[i][j]=0.0;

    /* loop over all nodes  */
    for (i=1; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        a = p->anc->index;   
        d = p->index;
        x = p->length;

        /*  start with distance from node to ancestor  */
        distsTemp[d][a] = x;
        distsTemp[a][d] = x;

        /* now revisit previously visited nodes, updating distances  */ 
        for (j=i-1; j>=0; j--)
            {
                k=(&t->nodes[j])->index;
                if (k==a) 
                { 
                    continue;
                }

                /*  dist from current node to prior node =  
                 *      dist from anc to prior node + dist from anc to current node  */
                distsTemp[d][k] = distsTemp[a][k] + x;
                distsTemp[k][d] = distsTemp[k][a] + x;

            }
        }

    /*  set the distances in the output array */
    for (i=0; i<(numExtNodes-1); i++) {
        for (j=i+1; j<numExtNodes; j++) {
            pd->dists[pairIdx(i,j,pd->nTaxa)]=distsTemp[i][j];
        }
    }

    /*  free temp array and return pointer to taxa pairwise distances */
    for (k=0;k<t->nNodes;k++)
        free(distsTemp[k]);
    free(distsTemp);

    return(NO_ERROR);
}


void CountFreqs(void) {

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

int CountDoublets(int nTaxa) {

    int s,i,j,k;
    int id1, id2;
 
    int I=4, J=4; 

    int nPairs = nTaxa*(nTaxa-1)/2;

    if (defMatrix == NO) 
    {
        MrBayesPrint("%s Matrix not Defined! \n", spacer);
    }

    if (defDoublets == YES) 
    {
        MrBayesPrint("%s Recomputing Doublets \n", spacer);
        free(pairwiseCounts);
    }

    pairwiseCounts = malloc(4*4*nPairs*sizeof(int));
    if (pairwiseCounts == NULL)
        {
        MrBayesPrint("%s Error Allocating Doublet Counts \n",spacer);
        return(ERROR);
        }

    // initialize doublet counts to 0
    for (k=0; k<nPairs; k++) 
        for (i=0; i<4; i++) 
             for (j=0; j<4; j++) {
                pairwiseCounts[tIdx(k,i,j,I,J)]=0;
             }

    for (i=1;i<nTaxa;i++) 
        {
        for (j=0;j<i;j++) 
            {
            k=pairIdx(i,j,nTaxa);

            for (s=0;s<numChar;s++) 
                {
                if (matrix[pos(i,s,numChar)]==GAP | matrix[pos(j,s,numChar)]==GAP)
                    continue;

                if (matrix[pos(i,s,numChar)]==MISSING | matrix[pos(j,s,numChar)]==MISSING)
                    continue;

                /* nucleotides at position x of sequences i & j   */        
                id2=toIdx(matrix[pos(j,s,numChar)]);
                id1=toIdx(matrix[pos(i,s,numChar)]);

                /* increment the count of that nucleotide pair, at pair k  */
                if (id1 == -1 || id2 == -1)
                    continue;

                pairwiseCounts[tIdx(k,id1,id2,I,J)]++;
                }
            }
        }

    defDoublets = YES;

    return(NO_ERROR);
}

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

void SetupQMat(MrBFlt **Q) {

    int i, j;

    if (defFreqs == NO) 
        {
        CountFreqs();
        }

    Q[0][1] = Q[1][0] = relRates[0];
    Q[0][2] = Q[2][0] = relRates[1];
    Q[0][3] = Q[3][0] = relRates[2];
    Q[1][2] = Q[2][1] = relRates[3];
    Q[1][3] = Q[3][1] = relRates[4];
    Q[2][3] = Q[3][2] = relRates[5];


    for (i=0;i<4;i++) {
        for (j=0;j<4;j++) {
            Q[i][j] = Q[i][j] * nucFreqs[j];
        }
    }

    // set diagonal elements: 
    double rowsum;
    for (i=0;i<4;i++) {
        rowsum=0.0;
        for (j=0;j<4;j++) {
            if (j!=i) rowsum+=Q[i][j];
        }
        Q[i][i]=-1*rowsum;
    }
} 

int DoPwSetParm (char *parmName, char *tkn)
{
    /*  char        *tempStr;  */
    MrBFlt      tempD;
    int j;

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Autoclose (autoClose) **********************************************************/
        if (!strcmp(parmName, "Dist"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                dist = tempD;
                MrBayesPrint ("%s   Setting distance for pw likelihood calcs to %lf\n", spacer, dist);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }        /* set Nowarnings (noWarn) **********************************************************/
        else if (!strcmp(parmName, "Alpha"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                alpha = tempD;
                MrBayesPrint ("%s   Setting alpha for pw likelihood calcs to %lf\n", spacer, alpha);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }        /* set Nowarnings (noWarn) **********************************************************/


        /* set Revmatpr (revMatPr) *********************************************************/
        else if (!strcmp(parmName, "Relrates"))
            {

            if (expecting == Expecting(EQUALSIGN)) 
                {
                expecting = Expecting(LEFTPAR);
                }

            else if (expecting == Expecting(LEFTPAR))
                {
                if (1) // TODO: add relrates to isargvalid (IsArgValid(tkn, tempStr) == NO_ERROR)

                    {
                    relRates[0] = relRates[1] = 1.0;
                    relRates[2] = relRates[3] = 1.0;
                    relRates[4] = relRates[5] = 1.0;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid Relrates argument\n", spacer);
                    return (ERROR);
                    }

                expecting  = Expecting(NUMBER);
                tempNumStates = 0;
                }

            else if (expecting == Expecting(NUMBER))
                {
                /* find out what type of prior is being set */
                /* find and store the number */
                sscanf (tkn, "%lf", &tempD);
                tempNum[tempNumStates++] = tempD;

                if (tempNumStates == 1)
                    expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                else if (tempNumStates < 6)
                    expecting  = Expecting(COMMA);
                else
                    expecting = Expecting(RIGHTPAR);
                }

            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }

            else if (expecting == Expecting(RIGHTPAR))
                {
                for (j=0; j<6; j++)
                    {
                    if (tempNumStates == 1)
                        relRates[j] = tempNum[0] / (MrBFlt) 6.0;
                    else
                        relRates[j] = tempNum[j];
                    }

                MrBayesPrint ("%s   Setting Relrates to (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
                    relRates[0], relRates[1], relRates[2],
                    relRates[3], relRates[4], relRates[5]);

                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);

                }
            else
                return (ERROR);

            } // closes parsing relrates param

        /* set Quitonerror (quitOnError) **************************************************/
        else
            {
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


MrBFlt* AllocateDists(int nDists) {

    MrBFlt* dists;
    dists=(MrBFlt*)malloc(nDists*sizeof(MrBFlt));
    return (dists);
}

int FreeDists(MrBFlt* dists) {

    free(dists);
    return (NO_ERROR);
}

int*** AllocateDoubletCounts(int nPairs) {

    int i;
    int ***doubletCounts=(int***)malloc(nPairs*sizeof(int**));
    for (i=0; i<nPairs; i++) {
        doubletCounts[i]=AllocateSquareIntegerMatrix(4);
    }
    return (doubletCounts);
}


void FreeDoubletCounts(int ***doubletCounts, int np) {
    int i;
    for (i=0; i<np; i++) {
        FreeSquareIntegerMatrix(doubletCounts[i]);
    }
    free(doubletCounts);
}


/* PrintNodes: Print a list of tree nodes, pointers and length */
void PrintPairwiseDists (PairwiseDists *pd)
{
    int i,j;
    int n;

    n = pd->nTaxa;
    /* printf ("tip1\ttip2\tdist\n"); */
    for (i=1; i<n; i++)
        {
        for (j=0; j<i; j++)
            {
            printf ("%d\t%d\t%f\n",
              i,j,
              pd->dists[pairIdx(i,j,n)]);
            }
        }

    printf ("\n");
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
                    doubP[dpIdx++] += bs[i] * tiP[index++]/((MrBFlt)m->numRateCats);
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
                    like=like*m->pwWeight;

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
            m->pwWeight = 2.0 / (numLocalTaxa * (numLocalTaxa - 1));
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
        MrBayesPrint("Calc pw weight for %d branch lengths \n", numBranches);

        V     = AllocateSquareDoubleMatrix(numBranches);
        Vinv  = AllocateSquareDoubleMatrix(numBranches);
             
        Vc    = AllocateSquareComplexMatrix(numBranches); 
        Vcinv = AllocateSquareComplexMatrix(numBranches);

        /*  * 
         *  Initialize necessary arrays:
         *  */
        /* initialize counts & index array */
        MrBayesPrint("Allocating counts and \n");
        MrBayesPrint("pairs %d \n", nPairs);
        MrBayesPrint("states %d \n", nStates);
        MrBayesPrint("splits %d \n", nSplits);

        countLen = (nSplits+1) * nPairs * nStates * nStates;
        MrBayesPrint("%s count array length: %d \n", spacer, countLen);

        counts = (int*) SafeMalloc( (nSplits+1) * nPairs * nStates * nStates * sizeof(int));
        if (counts == NULL)
            return(ERROR);

        MrBayesPrint("%s Done alloc counts \n", spacer);
        countIndex=(int**) SafeMalloc( (nSplits+1) * sizeof(int*));
        if (countIndex == NULL)
            return(ERROR);        

        for (i=0; i<(nSplits+1); i++)
            {
            countIndex[i]=(int*) SafeMalloc( nPairs * sizeof(int));
            if (countIndex[i] == NULL)
                return(ERROR);        
            }

        MrBayesPrint("%s Done alloc count index \n", spacer);

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
        MrBayesPrint("%s Done setting up count index \n", spacer);

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
        MrBayesPrint("done with inits \n");

        /*  *
         *  Done initializing
         *  */

        /* first calculate the n_ii for the division */
        MrBayesPrint("calc n_ii \n");
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
        MrBayesPrint("counts per data split \n");
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

//        MrBayesPrint("D1L\n");
//        for (i=0;i<4;i++) {
//            for (j=0;j<4;j++) {
//                MrBayesPrint(" %f ", D1L[i][j]);
//            }
//            MrBayesPrint("\n");
//        }
//
//        MrBayesPrint("D1LP\n");
//        for (i=0;i<4;i++) {
//            for (j=0;j<4;j++) {
//                MrBayesPrint(" %f ", D1LP[i][j]);
//            }
//            MrBayesPrint("\n");
//        }
        
//        MrBayesPrint(" %d ", n11[niiIndex[0][0]]);
//        MrBayesPrint(" %d \n", n10[niiIndex[0][0]]);
//
//        MrBayesPrint(" %d ", n11[niiIndex[4][0]]);
//        MrBayesPrint(" %d \n", n10[niiIndex[4][0]]);

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
            MrBayesPrint("Eigen %d = %f \n", i, eigvals[i]);
            eigsum += eigvals[i];
            eigsum2 += eigvals[i] * eigvals[i];
        }

        em = eigsum/(1.0*(numBranches)); 
        v = (eigsum * eigsum) / eigsum2;

        if (m->usePwWeights == 1)
            m->pwWeight=(1.0) / em;
        else if (m->usePwWeights == 2)  
            m->pwWeight=v / (1.0*numBranches+2.0*em);

        MrBayesPrint("%s pw weight: %f \n", spacer, m->pwWeight);

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
        if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!! \n");
           
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

    MrBayesPrint("PW Weights: GTR version \n \n \n ");
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
            m->pwWeight = 2.0 / (numLocalTaxa * (numLocalTaxa - 1));
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
        MrBayesPrint("Calc pw weight for %d branch lengths \n", numBranches);

        /*  * 
         *  Initialize necessary arrays:
         *  */
        /* initialize counts & index array */

        V     = AllocateSquareDoubleMatrix(numBranches);
        Vinv  = AllocateSquareDoubleMatrix(numBranches);
             
        Vc    = AllocateSquareComplexMatrix(numBranches); 
        Vcinv = AllocateSquareComplexMatrix(numBranches);

        MrBayesPrint("Allocating counts and \n");
        MrBayesPrint("pairs %d \n", nPairs);
        MrBayesPrint("states %d \n", nStates);
        MrBayesPrint("splits %d \n", nSplits);

        int countLen = (nSplits+1) * nPairs * nStates * nStates;
        MrBayesPrint("%s count array length: %d \n", spacer, countLen);

        counts = (int*) SafeMalloc( (nSplits+1) * nPairs * nStates * nStates * sizeof(int));
        if (counts == NULL)
            return(ERROR);

        MrBayesPrint("%s Done alloc counts \n", spacer);
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
                    counts[cI+dI] += 1;
                    counts[overallPwIdx+dI] += 1;
                    splitI+=1;
                    splitI = splitI%nSplits;

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

                            D1L[z][i]     += (counts[nidx+j] * t1);
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
        MrBFlt* dw= (MrBFlt *)SafeMalloc((size_t)numBranches*(sizeof(MrBFlt)));
        int*    iw= (int *)SafeMalloc((size_t)numBranches*(sizeof(int)));

        InvertMatrix(numBranches, H, dw,iw, Hinv);
        MultiplyMatrices(numBranches, Hinv, J, HiJ);

        //for (k=0; k<numBranches; k++)
        //    {
        //    MrBayesPrint("HiJ[%d][.] =",k);
        //    for (l=0; l<numBranches; l++)
        //        {
        //        MrBayesPrint("  % .3f", HiJ[k][l]);
        //        }
        //        MrBayesPrint("\n");
        //    }

        //MrBayesPrint("Is the problem here???? \n");
        int isComplex=GetEigens(numBranches,HiJ,eigvals,eigvalsc,V,Vinv,Vc,Vcinv);
        //MrBayesPrint("isComplex: %d \n", isComplex);

        MrBFlt eigsum=0.0 ;
        MrBFlt eigsum2=0.0 ;
        MrBFlt em=0.0;
        MrBFlt v=0.0;

        for (i=0; i<numBranches; i++) {
            //MrBayesPrint("Eigen %d = %f \n", i, eigvals[i]);
            eigsum += eigvals[i];
            eigsum2 += eigvals[i] * eigvals[i];
        }

        em = eigsum/(1.0*(numBranches)); 
        v = (eigsum * eigsum) / eigsum2;

        if (m->usePwWeights == 1)
            m->pwWeight=(1.0) / em;
        else if (m->usePwWeights == 2)  
            m->pwWeight=v / (1.0*numBranches+2.0*em);

        MrBayesPrint("%s pw weight: %f \n", spacer, m->pwWeight);

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

        /*  dists and probabilities */
        //MrBayesPrint("%s 1st free \n", spacer);        
        free(pwDists);
        //MrBayesPrint("%s 2nd free \n", spacer);
        //free(p_10 );
        //MrBayesPrint("%s 3nd free \n", spacer);
        //free(p_11 );
        //MrBayesPrint("%s 4th free \n", spacer);
        //free(p1_10); 
        //MrBayesPrint("%s 5th free \n", spacer);
        //free(p1_11);  
                      
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
