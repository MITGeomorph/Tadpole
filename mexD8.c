/*
 * =============================================================
 * mexD8.c - C MEX file to compute drainage area
 *
 * usage: [A, minima, flooded, W, Pav, Stack, Receiver, Zf, Lake] = mexD8(M,dy/dx,sinks,bvec,flood,P)
 *
 *    (A 7th input argument is accepted and ignored, so that mexD8, mexDinf
 *    and mexDms can be called with the same arguments.)
 *
 *    Input arguments:
 *      M:      a K x J matrix of elevations
 *      dy/dx:  ratio of grid spacings in the y and x directions
 *      sinks:  a K x J matrix of points where flow paths should terminate
 *      bvec:   a vector of codes for the [left, right, upper, lower]
 *              boundary conditions (0=fixed, 1=mirror, 2=periodic)
 *      flood:  a scalar. if 1, program will route flow through local 
 *              minima in the grid that do not drain to a sink. 
 *      P:      a K x J matrix of precipitation rates
 *
 *    Return arguments: 
 *      A:      a K x J matrix of total contributing areas in units of 
 *              cells.
 *      minima: a K x J matrix that is 1 where there are local minima or
 *              flat areas, zero otherwise.
 *      flooded:a matrix that is 1 where cells were originally in pits, zero otherwise.
 *      W:      a K x J x 8 array of flow routing weights
 *      Pav:    a K x J matrix of upstream-averaged precipitation rates
 *      Stack:  a (K*J) x 1 vector of 1-based (Matlab) linear indices of all
 *              cells, ordered so that every cell comes before all of its
 *              donors. Loop through it in reverse to visit every cell
 *              before any of its receivers.
 *      Receiver: a K x J x 8 array. Receiver(i,j,k) is the 1-based linear
 *              index of the cell that receives flow from (i,j) with weight
 *              W(i,j,k), or 0 where W(i,j,k) is 0.
 *      Zf:     a K x J matrix of elevations after Priority-Flood filling
 *              (equal to M where nothing was filled, or if flood = 0).
 *      Lake:   a K x J matrix. For each flooded cell, the 1-based linear
 *              index of the spill cell of its depression (all cells of one
 *              lake share it); 0 for cells that were not flooded.
 * 
 *    Important notes: 
 *      A must be multiplied by dx*dy to give the correct dimensions for
 *      the contributing area.
 *
 *      Fixed boundaries do not contribute flow.
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2006-2026 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>
#include "priorityflood.h"
#include "flowstack.h"

#define FIXED    0
#define MIRROR   1
#define PERIODIC 2
#define INF      1.0E10



// Calculate D8 (steepest descent) drainage directions
void D8Dir(const int K, const int J, double M[], double D[], double minima[], double minimaIdx[], int *numMin, double *a, double bdy[], double sinks[])
{

int i, j, k, Kj;
int di[8], di0[] = {0,-1,-1,-1,0,1,1,1};
int dj[8], dj0[] = {1,1,0,-1,-1,-1,0,1};
int skip[]={0,0,0,0,0,0,0,0};
int bl, br, bu, bd;
double invdx, invdy, invdiag, theS;
double e, s;
double invd[8];

invdx=1;
invdy=1/(*a);
invdiag=1/sqrt(1+(*a)*(*a));

invd[0] = invdx;
invd[1] = invdiag;
invd[2] = invdy;
invd[3] = invdiag;
invd[4] = invdx;
invd[5] = invdiag;
invd[6] = invdy;
invd[7] = invdiag;

// left
if (bdy[0] == 0) {
    bl=FIXED;
} else if (bdy[0] == 1) {
    bl=MIRROR;
} else if (bdy[0] == 2) {
    bl=PERIODIC;
}

// right
if (bdy[1] == 0) {
    br=FIXED;
} else if (bdy[1] == 1) {
    br=MIRROR;
} else if (bdy[1] == 2) {
    br=PERIODIC;
}

// upper
if (bdy[2] == 0) {
    bu=FIXED;
} else if (bdy[2] == 1) {
    bu=MIRROR;
} else if (bdy[2] == 2) {
    bu=PERIODIC;
}

// lower
if (bdy[3] == 0) {
    bd=FIXED;
} else if (bdy[3] == 1) {
    bd=MIRROR;
} else if (bdy[3] == 2) {
    bd=PERIODIC;
}



/* INTERIOR POINTS */

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}

for (j=1; j<(J-1); j++) { // Loop through columns
    
    Kj=K*j;
    
    for (i=1; i<(K-1); i++) { // Loop through rows
        
        if (!sinks[Kj+i]) { /* weights of sinks (which includes fixed boundaries) remain zero. Also, sinks and fixed boundaries are not defined as minima */
            e=M[Kj+i];
            theS=0;
            
            for (k=0; k<8; k++) {
                s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k]; // Could make this faster by only calculating slope if it's a downslope point
                if (s>theS) {
                    theS=s;
                    D[Kj+i]=k+1; // could make this faster by just recording k-steepest and only assigning D at the end 
                }
            }
            if (!theS) {
                // minima[Kj+i] = *numMin;
                // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
                // (*numMin)++;
                // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
                (*numMin)++;
                minima[Kj+i] = *numMin;
                minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
            }
        }
    }
}


/* BOUNDARIES, EXCLUDING CORNERS */

/* left */
j = 0;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (bl) {
case FIXED: // fixed left
case MIRROR: // mirror left
    skip[3]=1;
    skip[4]=1;
    skip[5]=1;
    break;
case PERIODIC: // periodic left
    dj[3]=J-1;
    dj[4]=J-1;
    dj[5]=J-1;
    break;
}

for (i=1; i<(K-1); i++) {
    if (!sinks[Kj+i]) { /* weights of sinks remain zero */
        e=M[Kj+i];
        theS=0;
        
        for (k=0; k<8; k++) {
            if (!skip[k]) {
                s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
                if (s>theS) {
                    theS=s;
                    D[Kj+i]=k+1;
                }
            }
        }
        if (!theS) {
            // minima[Kj+i] = *numMin;
            // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            // (*numMin)++;
            // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
            (*numMin)++;
            minima[Kj+i] = *numMin;
            minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
        } 
    }
}


/* right */
j = J-1;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (br) {
case FIXED: // fixed right
case MIRROR: // mirror right
    skip[7]=1;
    skip[0]=1;
    skip[1]=1;
    break;
case PERIODIC: // periodic right
    dj[7]=1-J;
    dj[0]=1-J;
    dj[1]=1-J;
    break;
}

for (i=1; i<(K-1); i++) { 
    if (!sinks[Kj+i]) { /* weights of sinks remain zero */
        e=M[Kj+i];
        theS=0;

        for (k=0; k<8; k++) {
            if (!skip[k]) {
                s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
                if (s>theS) {
                    theS=s;
                    D[Kj+i]=k+1;
                }
            }
        }
        if (!theS) {
            // minima[Kj+i] = *numMin;
            // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            // (*numMin)++;
            // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
            (*numMin)++;
            minima[Kj+i] = *numMin;
            minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
        } 
    }
}


/* upper */
i = 0;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (bu) {
case FIXED: // fixed upper
case MIRROR: // mirror upper
    skip[1]=1;
    skip[2]=1;
    skip[3]=1;
    break;
case PERIODIC: // periodic upper
    di[1]=K-1;
    di[2]=K-1;
    di[3]=K-1;
    break;
}

for (j=1; j<(J-1); j++) {
    Kj = K*j;
    if (!sinks[Kj+i]) { /* weights of sinks remain zero */
        e=M[Kj+i];
        theS=0;

        for (k=0; k<8; k++) {
            if (!skip[k]) {
                s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
                if (s>theS) {
                    theS=s;
                    D[Kj+i]=k+1;
                }
            }
        }
        if (!theS) {
            // minima[Kj+i] = *numMin;
            // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            // (*numMin)++;
            // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
            (*numMin)++;
            minima[Kj+i] = *numMin;
            minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
        }         
    }
}


/* lower */
i = K-1;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (bd) {
case FIXED: // fixed lower
case MIRROR: // mirror lower
    skip[5]=1;
    skip[6]=1;
    skip[7]=1;
    break;
case PERIODIC: // periodic lower
    di[5]=1-K;
    di[6]=1-K;
    di[7]=1-K;
    break;
}

for (j=1; j<(J-1); j++) {
    Kj = K*j;
    if (!sinks[Kj+i]) { /* weights of sinks remain zero */
        e=M[Kj+i];
        theS=0;

        for (k=0; k<8; k++) {
            if (!skip[k]) {
                s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
                if (s>theS) {
                    theS=s;
                    D[Kj+i]=k+1;
                }
            }
        }
        if (!theS) {
            // minima[Kj+i] = *numMin;
            // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            // (*numMin)++;
            // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
            (*numMin)++;
            minima[Kj+i] = *numMin;
            minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
        }         
    }
}


/* CORNERS */

/* UL */
i=0; 
j=0;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (bl) {
case FIXED: // fixed left
case MIRROR: // mirror left
    skip[3]=1;
    skip[4]=1;
    skip[5]=1;
    break;
case PERIODIC: // periodic left
    dj[3]=J-1;
    dj[4]=J-1;
    dj[5]=J-1;
    break;
}

switch (bu) {
case FIXED: // fixed upper
case MIRROR: // mirror upper
    skip[1]=1;
    skip[2]=1;
    skip[3]=1;
    break;
case PERIODIC: // periodic upper
    di[1]=K-1;
    di[2]=K-1;
    di[3]=K-1;
    break;
}

if (!sinks[Kj+i]) { /* weights of sinks remain zero */
    e=M[Kj+i];
    theS=0;

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
            if (s>theS) {
                theS=s;
                D[Kj+i]=k+1;
            }
        }
    }
    if (!theS) {
        // minima[Kj+i] = *numMin;
        // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        // (*numMin)++;
        // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
        (*numMin)++;
        minima[Kj+i] = *numMin;
        minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
    } 
}


/* UR */
i=0;
j=J-1;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (br) {
case FIXED: // fixed right
case MIRROR: // mirror right
    skip[7]=1;
    skip[0]=1;
    skip[1]=1;
    break;
case PERIODIC: // periodic right
    dj[7]=1-J;
    dj[0]=1-J;
    dj[1]=1-J;
    break;
}

switch (bu) {
case FIXED: // fixed upper
case MIRROR: // mirror upper
    skip[1]=1;
    skip[2]=1;
    skip[3]=1;
    break;
case PERIODIC: // periodic upper
    di[1]=K-1;
    di[2]=K-1;
    di[3]=K-1;
    break;
}

if (!sinks[Kj+i]) { /* weights of sinks remain zero */
    e=M[Kj+i];
    theS=0;

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
            if (s>theS) {
                theS=s;
                D[Kj+i]=k+1;
            }
        }
    }
    if (!theS) {
        // minima[Kj+i] = *numMin;
        // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        // (*numMin)++;
        // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
        (*numMin)++;
        minima[Kj+i] = *numMin;
        minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
    } 
}


/* LL */
i=K-1;
j=0;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (bl) {
case FIXED: // fixed left
case MIRROR: // mirror left
    skip[3]=1;
    skip[4]=1;
    skip[5]=1;
    break;
case PERIODIC: // periodic left
    dj[3]=J-1;
    dj[4]=J-1;
    dj[5]=J-1;
    break;
}

switch (bd) {
case FIXED: // fixed lower
case MIRROR: // mirror lower
    skip[5]=1;
    skip[6]=1;
    skip[7]=1;
    break;
case PERIODIC: // periodic lower
    di[5]=1-K;
    di[6]=1-K;
    di[7]=1-K;
    break;
}

if (!sinks[Kj+i]) { /* weights of sinks remain zero */
    e=M[Kj+i];
    theS=0;

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
            if (s>theS) {
                theS=s;
                D[Kj+i]=k+1;
            }
        }
    }
    if (!theS) {
        // minima[Kj+i] = *numMin;
        // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        // (*numMin)++;
        // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
        (*numMin)++;
        minima[Kj+i] = *numMin;
        minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
    } 
}


/* LR */
i=K-1;
j=J-1;
Kj=K*j;

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}
for (k=0; k<8; k++) {skip[k]=0;}

switch (br) {
case FIXED: // fixed right
case MIRROR: // mirror right
    skip[7]=1;
    skip[0]=1;
    skip[1]=1;
    break;
case PERIODIC: // periodic right
    dj[7]=1-J;
    dj[0]=1-J;
    dj[1]=1-J;
    break;
}

switch (bd) {
case FIXED: // fixed lower
case MIRROR: // mirror lower
    skip[5]=1;
    skip[6]=1;
    skip[7]=1;
    break;
case PERIODIC: // periodic lower
    di[5]=1-K;
    di[6]=1-K;
    di[7]=1-K;
    break;
}

if (!sinks[Kj+i]) { /* weights of sinks remain zero */
    e=M[Kj+i];
    theS=0;

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            s = (e - M[K*(j+dj[k])+(i+di[k])])*invd[k];
            if (s>theS) {
                theS=s;
                D[Kj+i]=k+1;
            }
        }
    }
    if (!theS) {
        // minima[Kj+i] = *numMin;
        // minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        // (*numMin)++;
        // The above code contains a bug: numMin starts at zero, so for the first minimum found, minima[Kj+i] = 0, which will make it indistinguishable from the non-minimum elements in minima[]. The code below fixes the bug.
        (*numMin)++;
        minima[Kj+i] = *numMin;
        minimaIdx[(*numMin)-1] = Kj+i; // Add this location to list of indices of minima
    } 
}


} /* end D8Dir function */


void FixMirrorWeights(const int K, const int J, double W[], double bdy[]) {
    
/* Halve the weights of flow from a cell on a mirror boundary to a cell off that
   boundary. Corner cells are included: a corner between two mirror boundaries is
   halved once for each of those boundaries its receiver is not on (so a diagonal
   weight into the interior is quartered), and a corner between a mirror boundary
   and another boundary type is halved like the rest of the mirror boundary. */
int i, j, k1, k2, k3;

/* left */
if (bdy[0] == 1) {
    j = 0;
    k1 = 7;
    k2 = 0;
    k3 = 1;
    
    for (i=0; i<K; i++) {
        W[k1*K*J+K*j+i] *= 0.5;
        W[k2*K*J+K*j+i] *= 0.5;
        W[k3*K*J+K*j+i] *= 0.5;
    }
}

/* right */
if (bdy[1] == 1) {
    j = J-1;
    k1 = 3;
    k2 = 4;
    k3 = 5;
    
    for (i=0; i<K; i++) {
        W[k1*K*J+K*j+i] *= 0.5;
        W[k2*K*J+K*j+i] *= 0.5;
        W[k3*K*J+K*j+i] *= 0.5;
    }
}


/* upper */
if (bdy[2] == 1) {
    i = 0;
    k1 = 5;
    k2 = 6;
    k3 = 7;
    
    for (j=0; j<J; j++) {
        W[k1*K*J+K*j+i] *= 0.5;
        W[k2*K*J+K*j+i] *= 0.5;
        W[k3*K*J+K*j+i] *= 0.5;
    }
}


/* lower */
if (bdy[3] == 1) {
    i = K-1;
    k1 = 1;
    k2 = 2;
    k3 = 3;
    
    for (j=0; j<J; j++) {
        W[k1*K*J+K*j+i] *= 0.5;
        W[k2*K*J+K*j+i] *= 0.5;
        W[k3*K*J+K*j+i] *= 0.5;
    }
}
    
} /* end FixMirrorWeights() */


void ReplaceWeights(const int K, const int J, double W[], double D[]) {

    int i, j, k;

    for (i=0; i<K; i++) {
        for (j=0; j<J; j++) {
            if ( (k = (int)D[K*j+i]) ) { // if i,j has a nonzero D8 direction
                W[(k-1)*K*J+K*j+i]=1; // copy the D8 weight
            }
        }
    }

} // end ReplaceWeights()



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        
    mxArray *Dptr, *minimaIdxptr;
    double *Zf = NULL, *Lake = NULL;
    double *M, *Mfl, *P, *Pav, *Stack = NULL, *Receiver = NULL, *minimaIdx, *D, *F, *W, *A, *minima, *a, *fl, *b, *s;
    int n, K, J, *numMin, ndims=3;
    mwSize dims[]={0, 0, 8};
    
    if (nrhs < 6 || nrhs > 7)
        mexErrMsgIdAndTxt("mexD8:nrhs", "Usage: [A, minima, flooded, W, Pav, Stack, Receiver, Zf, Lake] = mexD8(M,dy/dx,sinks,bvec,flood,P[,p])");
    if (nlhs > 9)
        mexErrMsgIdAndTxt("mexD8:nlhs", "mexD8 returns at most 9 outputs.");

    // Get pointers to inputs
    M = (double *)mxGetPr(prhs[0]); // elevations
    a = (double *)mxGetPr(prhs[1]); // dy/dx
    s = (double *)mxGetPr(prhs[2]); /* sinks where recursive upslope area function should be called */
    b = (double *)mxGetPr(prhs[3]); /* [bleft bright bupper blower] */
    fl = (double *)mxGetPr(prhs[4]); // do or don't do flooding
    P = (double *)mxGetPr(prhs[5]); // precipitation rates
    
    // Get dimensions of input matrix of elevations
    K = (int)mxGetM(prhs[0]);
    J = (int)mxGetN(prhs[0]);
    dims[0] = K;
    dims[1] = J;    
    
    // Create arrays for return arguments
    minima = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // 1 if an element is a local minimum, zero otherwise
    F = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix that indicates what was flooded
    if (nlhs > 7) Zf = (double *)mxGetPr(plhs[7]= mxCreateDoubleMatrix(K, J, mxREAL)); // filled surface
    if (nlhs > 8) Lake = (double *)mxGetPr(plhs[8]= mxCreateDoubleMatrix(K, J, mxREAL)); // lake labels (spill cell index), 0 if not flooded

    minimaIdx = (double *)mxGetPr(minimaIdxptr= mxCreateDoubleMatrix(K, J, mxREAL)); 
    D = (double *)mxGetPr(Dptr= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of D8 drainage directions
    
    numMin = mxCalloc(1, sizeof(int)); // allocate memory for this 1-element integer array
    *numMin = 0;
    
    
    /* FLOODING */
    if (*fl) {  // If user requested flooding
        
        // Fill depressions and flats with Priority-Flood+epsilon, so every cell other than sinks and fixed boundaries has a strictly lower neighbor
        Mfl = (double *)mxMalloc(K*J*sizeof(double)); // the post-flood elevations
        PriorityFloodEpsilon(K, J, M, Mfl, s, NULL, b, PF_EPS_AUTO, NULL, Lake);
        if (Zf) memcpy(Zf, Mfl, K*J*sizeof(double)); // return the filled surface

        // Flag the cells that were raised by flooding
        for (n=0; n<(K*J); n++) {
            F[n] = (Mfl[n] > M[n]);
        }

        // Get D8 drainage directions using the flooded elevations
        D8Dir(K, J, Mfl, D, minima, minimaIdx, numMin, a, b, s);

        /* free memory */
        mxFree(Mfl);

    } else { // flooding was not requested
                
        // Get D8 drainage directions without flooding
        D8Dir(K, J, M, D, minima, minimaIdx, numMin, a, b, s);
        
        if (Zf) memcpy(Zf, M, K*J*sizeof(double)); // nothing was filled
        
        // No cells were flooded
        for (n = 0; n < K*J; n++) {
            F[n] = 0;
        }
                
    }
    
    // The weights array
    W = (double *)mxGetPr(plhs[3]= mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL));

    // copy D8 directions into the weights array
    ReplaceWeights(K, J, W, D);
    
    // Free memory
    mxFree(numMin);
    mxDestroyArray(Dptr);
    mxDestroyArray(minimaIdxptr);
    
    /* If necessary, correct weights for points on mirror boundaries */
    /* (A point on a mirrored boundary donating flow to a point off the boundary should only donate half its area) */
    if (b[0]==1 || b[1]==1 || b[2]==1 || b[3]==1) {
        FixMirrorWeights(K, J, W, b);
    }
    
    
    // Create arrays for additional return arguments 
    A = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of TCAs
    Pav = (double *)mxGetPr(plhs[4]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of upstream-averaged precipitation rates
    if (nlhs > 5) Stack = (double *)mxGetPr(plhs[5]= mxCreateDoubleMatrix(K*J, 1, mxREAL)); // 1-based linear indices in stack order
    if (nlhs > 6) Receiver = (double *)mxGetPr(plhs[6]= mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL)); // 1-based linear indices of receivers
    
    // Get contributing area for each element, and the stack and receivers if requested
    FlowAccumulate(K, J, W, b, P, A, Pav, Stack, Receiver);

    // Calculate upstream-averaged precipitation rates by dividing the weighted sum of precip rates by the sum of the weights
    for (n=0; n<(K*J); n++) {
        Pav[n] = Pav[n]/A[n];
        // Pav[n] = Mfl[n]; // debugging: copy the post-flood elevations into Pav!
    }

    
} // End mexFunction