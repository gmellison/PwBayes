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

#include "bayes.h"
#include "pairwise.h"
#include "utils.h"
#include "command.h"

#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

// global variables for pairwise likelihoods:
// global variables for pw likelihoods:
MrBFlt nucFreqs[4];
int    doubletCounts[4][4];
MrBFlt doubletFreqs[4][4];
int    tripleCounts[64];
MrBFlt relRates[6];
int tempNumStates;

int defFreqs = NO;
int defDoublets = NO;
int defTriples = NO;

MrBFlt alpha=0.1;
MrBFlt dist=0.1;

// local prototypea:
int   toIdx(int x);
int   DoEstQPairwise(void);

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


void CountDoublets(void) {

    int s,i,j;
    int id1, id2;
    
    // initialize doublet counts to 0
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++) 
                doubletCounts[i][j]=0;

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
        for (j=0;j<4;j++) {
            doubletFreqs[i][j] = (doubletCounts[i][j] + doubletCounts[j][i]) / (2.0 * numTaxa * (numTaxa-1)*numChar);
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


void CountTriples(void) {

    int s, seqi, seqj, seqk;
    int i,j,k;

    // initialize triplet counts to 0
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++) 
                for (k=0;k<4;k++)
                    tripleCounts[triplePos(i,j,k)]=0;

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

int DoEstQPairwise(void) {

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }
           
    if (defFreqs == NO) 
    {
        CountFreqs(); // sets nucleotide frequencies in global object 'nucFreqs'
    }

    if (defDoublets == NO) 
    {
        CountDoublets(); // sets nucleotide frequencies in global object 'doubletFreqs'
    }

    MrBayesPrint("Running 'Estqpairwise' command \n");

    int i,j;//,s;//,id2,id2;
    double nu=0.0;

    // use mb machinery to compute eigens
    // set up matrices 
    MrBFlt **V         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Vinv      = AllocateSquareDoubleMatrix(4);
    MrBFlt **LaLog     = AllocateSquareDoubleMatrix(4);

    MrBComplex **Vc    = AllocateSquareComplexMatrix(4); 
    MrBComplex **Vcinv = AllocateSquareComplexMatrix(4);
    MrBFlt la[4]; 
    MrBFlt laC[4];

    MrBFlt **F         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Q         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Temp      = AllocateSquareDoubleMatrix(4);

    // set up matrix of empirical transitions 
    //   probabilities, 'Q' (will be modified in place)  
    for (i=0;i<4;i++) {
        for (j=0;j<4;j++){
            F[i][j] = doubletFreqs[i][j] / nucFreqs[j];
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

    // calculate the matrix log log(e^Qt) = V log[La] V^-1) 
    MultiplyMatrices(4,V,LaLog,Temp); 
    MultiplyMatrices(4,Temp,Vinv,Q); // Q is the ptr to resulting matrix
    FreeSquareDoubleMatrix(Temp);
    FreeSquareDoubleMatrix(F);

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
        nu += -1.0 * Q[i][i] * nucFreqs[i];

    MultiplyMatrixByScalar(4, Q, 1.0/nu, Q);

    MrBayesPrint("%s Tau Hat: %f \n", spacer, nu);

    // get the individual rates: 
    double rates[6]; 
    rates[0] = Q[0][1];
    rates[1] = Q[0][2];
    rates[2] = Q[0][3];
    rates[3] = Q[1][2];
    rates[4] = Q[1][3];
    rates[5] = Q[2][3];

    MrBayesPrint("%s Estimated Relative Rates:", spacer);

    for (int r=0;r<6;r++)
            MrBayesPrint(" %f", rates[r]);
    MrBayesPrint("\n");

    // free all matrices 
    FreeSquareDoubleMatrix(V);
    FreeSquareDoubleMatrix(Vinv);
    FreeSquareDoubleMatrix(Q);
    FreeSquareComplexMatrix(Vc);
    FreeSquareComplexMatrix(Vcinv);

    return(NO_ERROR);
}

void SetupQMat(MrBFlt **Q, MrBFlt *rates, MrBFlt *freqs) {

    int i, j;

    // set up Q matrix, using input rates & empirical freqs:
    Q[0][1] = Q[1][0] = rates[0];
    Q[0][2] = Q[2][0] = rates[1];
    Q[0][3] = Q[3][0] = rates[2];
    Q[1][2] = Q[2][1] = rates[3];
    Q[1][3] = Q[3][1] = rates[4];
    Q[2][3] = Q[3][2] = rates[5];


    for (i=0;i<4;i++) {
        for (j=0;j<4;j++) {
            Q[i][j] = Q[i][j] * freqs[j];
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


int DoPairwiseLogLike(void) {


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
        CountDoublets();
    }

    MrBayesPrint("Running 'Pairwiseloglike' command.\n");

    int i,j;
    int ri;
    MrBFlt tp;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

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

    SetupQMat(Q,relRates,nucFreqs);

    //MrBFlt probsByRateCat[4];
    MrBFlt ll=0.0;

    // calculate the site pattern probability for each 
    //   dg4 rate category
    for (ri=0; ri<4; ri++) {
    
        // probability transition matrix for site rate i:
        MultiplyMatrixByScalar(4, Q, dist*dg4[ri], Qtausr);  

        int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
        if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!!");
           
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
            ll += doubletCounts[i][j] * log(tp);
            
        }
    }

    FreeSquareDoubleMatrix(TempMat);
    for (int ri=0;ri<4;ri++)
        FreeSquareDoubleMatrix(pTrans[ri]);

    MrBayesPrint("Pairwise Log-Likelihood: %f \n", ll);
    return(NO_ERROR);

}


int DoTripletLogLike(void) {

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }

    if (defFreqs==NO) 
    {
        CountFreqs();
    }

    if (defTriples==NO)
    {
        CountTriples();
    }

    // temp: just set distance, alpha:
    double bl=0.1;
    double al=0.5;
   
    MrBFlt dg4[4];
    DiscreteGamma(dg4,al,al,4,1);

    MrBFlt pii[4];
    MrBFlt pij[4];

    for (int i=0; i<4; i++) {
        MrBFlt etr = exp(-(4.0/3) * bl * dg4[i]);
        pii[i]=(1.0/4) + (3.0/4) * etr;
        pij[i]=(1.0/4) - (1.0/4) * etr;
    }

    MrBFlt tripleProbs[64], transitionProbs[64];

    int i,j,k, ri;

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
                // no need to calculate probs for not present in data
                if (tripleCounts[triplePos(i,j,k)] == 0) continue;  

                // ensure triple site pattern prob is init to 0:
                tripleProbs[triplePos(i,j,k)] = 0.0;

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
                    tripleProbs[triplePos(i,j,k)]+=0.25*probsByRateCat[ri];
                }

                ll += tripleCounts[triplePos(i,j,k)] * log(tripleProbs[triplePos(i,j,k)]);
            }
        }
    }

    MrBayesPrint("Triplet Quasi Log-Likelihood for Data: %f", ll);
    return(NO_ERROR);

}


int DoPwSetParm (char *parmName, char *tkn)
{
    char        *tempStr;
    MrBFlt      tempD;
    int flag;
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
                    flag = 0;
                    relRates[0] = relRates[1] = 1.0;
                    relRates[2] = relRates[3] = 1.0;
                    relRates[4] = relRates[5] = 1.0;
                    flag = 1;
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


/* PrintNodes: Print a list of tree nodes, pointers and length */
void PrintPairwiseDists (PolyTree *t)
{
    int         i,j,a,d,k;
    PolyNode    *p,*e1,*e2;
    double      x;

    /* printf ("tip1\ttip2\tdist\n"); */

    MrBFlt **dists = AllocateSquareDoubleMatrix(t->nNodes);

    for (i=0; i<t->nNodes; i++) 
        for (j=0; j<t->nNodes; j++)
                dists[i][j]=0.0;

    for (i=1; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        a = p->anc->index;   
        d = p->index;
        x = p->length;

        dists[a][d] = dists[a][d] += x;
        
        for (j=i-1; j>=0; j--)
            {
                k=(&t->nodes[j])->index;
                if (k==a) continue;
                dists[d][k] = dists[k][d] += x;
            }
        }

    for (i=1; i<t->nNodes-1; i++)
        {
        for (j=0; j<i; j++)
            {
            printf ("%d\t%d\t%f\n",
            i,j,
            dists[i][j]);
            }
        }

    FreeSquareDoubleMatrix(dists);
    printf ("\n");
}


