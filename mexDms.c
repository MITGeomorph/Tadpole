 /*
 * =============================================================
 * mexDms.c - C MEX file to compute drainage area
 *
 * usage: [A, minima, flooded, W] = mexDms(M,dy/dx,sinks,bvec,flood)
 *
 *    Input arguments:
 *      M:      a K x J matrix of elevations
 *      dy/dx:  ratio of grid spacings in the y and x directions
 *      sinks:  a K x J matrix of points where flow paths should terminate
 *      bvec:   a vector of codes for the [left, right, upper, lower]
 *              boundary conditions (0=fixed, 1=mirror, 2=periodic)
 *      flood:  a scalar. if 1, program will route flow through local 
 *              minima in the grid that do not drain to a sink. 
 *
 *    Return arguments: 
 *      A:      a K x J matrix of total contributing areas in units of 
 *              cells, calculated using the multi-slope algorithm.
 *      minima: a K x J matrix that is 1 where there are local minima or
 *              flat areas, zero otherwise.
 *      flooded:a matrix that is 1 where M was flooded, zero otherwise.
 *      W:      a K x J x 8 array of flow routing weights
 *
 *    Important notes: 
 *      A must be multiplied by dx*dy to give the correct dimensions for
 *      the contributing area.
 *
 *      Fixed boundaries do not contribute flow.
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2006-2011 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

#define FIXED    0
#define MIRROR   1
#define PERIODIC 2
#define INF      1.0E10
#define MAX_NEIGHBOR_BASINS 1000

typedef struct {
    int nnb, lpp, inb[MAX_NEIGHBOR_BASINS], ipp[MAX_NEIGHBOR_BASINS]; // number of neighboring basins, list position of the lowest pp in this structure, index of neighboring basin, index of pp
    double zpp[MAX_NEIGHBOR_BASINS]; // elevation of pp
} pourpoint;


// Calculate D8 (steepest descent) drainage directions
void D8Dir(const int K, const int J, double M[], double D[], double minima[], double minimaIdx[], int *numMin, double *a, double bdy[], double sinks[])
{

int i, j, k, Kj, theD;
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
                minima[Kj+i] = *numMin;
                minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
                (*numMin)++;
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
            minima[Kj+i] = *numMin;
            minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            (*numMin)++;
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
            minima[Kj+i] = *numMin;
            minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            (*numMin)++;
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
            minima[Kj+i] = *numMin;
            minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            (*numMin)++;
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
            minima[Kj+i] = *numMin;
            minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
            (*numMin)++;
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
        minima[Kj+i] = *numMin;
        minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        (*numMin)++;
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
        minima[Kj+i] = *numMin;
        minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        (*numMin)++;
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
        minima[Kj+i] = *numMin;
        minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        (*numMin)++;
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
        minima[Kj+i] = *numMin;
        minimaIdx[*numMin] = Kj+i; // Add this location to list of indices of minima
        (*numMin)++;
    } 
}


} /* end D8Dir function */



void GetBasin(const int i, const int j, const int label, const int K, const int J, double D[], double Basin[], double bdy[])
{

int k, p, q, bl, br, bu, bd;
int r[] = {5,6,7,8,1,2,3,4};
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};

if (Basin[K*j+i]==label) {return;} // Bail out if this cell is already part of the basin

Basin[K*j+i]=label; // Now it is part of the basin

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


if (i==0) {
    
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

} else if (i==K-1) {

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


}

if (j==0) {
    
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
    
} else if (j==J-1) {

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

}


/* loop through neighbors */
for (k=0; k<8; k++) {

    if (!skip[k]) {
    
        p = i + di[k];
        q = j + dj[k];

        if (D[K*q+p] == r[k]) { // if neighbor k drains to current cell
            GetBasin(p,q,label,K,J,D,Basin,bdy);
        }

    } 

}

} // End GetBasin



void Repaint(const int i, const int j, const int nold, const int nnew, const int K, const int J, double Basin[], double bdy[])
{

int k, p, q, bl, br, bu, bd;
int r[] = {5,6,7,8,1,2,3,4};
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};

Basin[K*j+i] = nnew; // Now it is part of the new basin

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


if (i==0) {
    
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

} else if (i==K-1) {

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


}

if (j==0) {
    
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
    
} else if (j==J-1) {

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

}


/* loop through neighbors */
for (k=0; k<8; k++) {

    if (!skip[k]) {
    
        p = i + di[k];
        q = j + dj[k];

        if (Basin[K*q+p] == nold) { // if neighbor k is part of the old basin
            Repaint(p,q,nold,nnew,K,J,Basin,bdy);
        }

    } 

}


} // End Repaint


void MemCheck(pourpoint pp[], int n) {
    if (pp[n].nnb > MAX_NEIGHBOR_BASINS) {
        mexErrMsgTxt("Memory overrun in mexDms. Try increasing MAX_NEIGHBOR_BASINS and recompiling, or using a smaller grid with fewer local minima.");
    }
} // End MemCheck


void GetPourPoints(pourpoint pp[],const int K, const int J, double M[], double Basin[], double minima[], double minimaIdx[], int *numMin, double sinks[], double bdy[])
{

int i, j, k, Kj, p, q, c, bl, br, bu, bd;
int bas, nbas, ncand, midx, nmidx, lowestppyet, ppexists, pploc;
int di[8], di0[] = {0,-1,-1,-1,0,1,1,1};
int dj[8], dj0[] = {1,1,0,-1,-1,-1,0,1};
int skip[]={0,0,0,0,0,0,0,0};

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
        
        bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to
        
        if (!sinks[bas]) { // only process points in closed basins
            
            for (k=0; k<8; k++) { // loop through neighbors
                p = i + di[k];
                q = j + dj[k];
                
                nbas = (int)Basin[K*q+p];
                
                if (nbas != bas) { // if the neighbor is part of a different basin
                    if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                        ncand = K*q+p; // the neighbor is a candidate pour point
                    } else { // the neighbor has equal or lower elevation
                        ncand = K*j+i; // the present point is a candidate pour point
                    }
                    
                    // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                    midx = (int)minima[bas]; // the list position of the minimum that this basin drains to
                    
                    lowestppyet = 0;
                    ppexists = 0;
                    
                    for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                        if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                            ppexists = 1;
                            pploc = c;
                            if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                                lowestppyet = 1;
                            }
                            break;
                        }
                    }
                    
                    if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                        pp[midx].inb[pp[midx].nnb] = nbas;
                        pp[midx].ipp[pp[midx].nnb] = ncand;
                        pp[midx].zpp[pp[midx].nnb] = M[ncand];
                        (pp[midx].nnb)++; MemCheck(pp,midx);
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            pp[nmidx].inb[pp[nmidx].nnb] = bas;
                            pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                            pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                            (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                        }
                    } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                        pp[midx].inb[pploc] = nbas;
                        pp[midx].ipp[pploc] = ncand;
                        pp[midx].zpp[pploc] = M[ncand];
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            // search for the position of the existing PP in the neighbor basin's list
                            c = 0;
                            while (pp[nmidx].inb[c] != bas) { c++; }
                            pp[nmidx].inb[c] = bas;
                            pp[nmidx].ipp[c] = ncand;
                            pp[nmidx].zpp[c] = M[ncand];
                        }
                    } // otherwise, there was an existing pour point that was lower than the candidate; do nothing
                }
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
    
    bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to
    
    if (!sinks[bas]) { /* weights of sinks remain zero */
        
        for (k=0; k<8; k++) {
            if (!skip[k]) {
                p = i + di[k];
                q = j + dj[k];
                
                nbas = (int)Basin[K*q+p];
                
                if (nbas != bas) { // if the neighbor is part of a different basin
                    if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                        ncand = K*q+p; // the neighbor is a candidate pour point
                    } else { // the neighbor has equal or lower elevation
                        ncand = K*j+i; // the present point is a candidate pour point
                    }
                    
                    // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                    midx = (int)minima[bas]; // the list position of the minimum that this basin drains to
                    
                    lowestppyet = 0;
                    ppexists = 0;
                    for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                        if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                            ppexists = 1;
                            pploc = c;
                            if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                                lowestppyet = 1;
                            }
                            break;
                        }
                    }
                    
                    if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                        pp[midx].inb[pp[midx].nnb] = nbas;
                        pp[midx].ipp[pp[midx].nnb] = ncand;
                        pp[midx].zpp[pp[midx].nnb] = M[ncand];
                        (pp[midx].nnb)++; MemCheck(pp,midx);
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            pp[nmidx].inb[pp[nmidx].nnb] = bas;
                            pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                            pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                            (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                        }
                    } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                        pp[midx].inb[pploc] = nbas;
                        pp[midx].ipp[pploc] = ncand;
                        pp[midx].zpp[pploc] = M[ncand];
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            // search for the position of the existing PP in the neighbor basin's list
                            c = 0;
                            while (pp[nmidx].inb[c] != bas) { c++; }
                            pp[nmidx].inb[c] = bas;
                            pp[nmidx].ipp[c] = ncand;
                            pp[nmidx].zpp[c] = M[ncand];
                        }
                    }
                }
            }
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
    
    bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to
    
    if (!sinks[bas]) { /* weights of sinks remain zero */
        
        for (k=0; k<8; k++) {
            if (!skip[k]) {
                p = i + di[k];
                q = j + dj[k];
                
                nbas = (int)Basin[K*q+p];
                
                if (nbas != bas) { // if the neighbor is part of a different basin
                    if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                        ncand = K*q+p; // the neighbor is a candidate pour point
                    } else { // the neighbor has equal or lower elevation
                        ncand = K*j+i; // the present point is a candidate pour point
                    }
                    
                    // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                    midx = (int)minima[bas]; // the list position of the minimum that this basin drains to
                    
                    lowestppyet = 0;
                    ppexists = 0;
                    for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                        if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                            ppexists = 1;
                            pploc = c;
                            if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                                lowestppyet = 1;
                            }
                            break;
                        }
                    }
                    
                    if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                        pp[midx].inb[pp[midx].nnb] = nbas;
                        pp[midx].ipp[pp[midx].nnb] = ncand;
                        pp[midx].zpp[pp[midx].nnb] = M[ncand];
                        (pp[midx].nnb)++; MemCheck(pp,midx);
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            pp[nmidx].inb[pp[nmidx].nnb] = bas;
                            pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                            pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                            (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                        }
                    } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                        pp[midx].inb[pploc] = nbas;
                        pp[midx].ipp[pploc] = ncand;
                        pp[midx].zpp[pploc] = M[ncand];
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            // search for the position of the existing PP in the neighbor basin's list
                            c = 0;
                            while (pp[nmidx].inb[c] != bas) { c++; }
                            pp[nmidx].inb[c] = bas;
                            pp[nmidx].ipp[c] = ncand;
                            pp[nmidx].zpp[c] = M[ncand];
                        }
                    }
                }
            }
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
    
    bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to
    
    if (!sinks[bas]) { /* weights of sinks remain zero */
        
        for (k=0; k<8; k++) {
            if (!skip[k]) {
                p = i + di[k];
                q = j + dj[k];
                
                nbas = (int)Basin[K*q+p];
                
                if (nbas != bas) { // if the neighbor is part of a different basin
                    if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                        ncand = K*q+p; // the neighbor is a candidate pour point
                    } else { // the neighbor has equal or lower elevation
                        ncand = K*j+i; // the present point is a candidate pour point
                    }
                    
                    // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                    midx = (int)minima[bas]; // the list position of the minimum that this basin drains to
                    
                    lowestppyet = 0;
                    ppexists = 0;
                    for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                        if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                            ppexists = 1;
                            pploc = c;
                            if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                                lowestppyet = 1;
                            }
                            break;
                        }
                    }
                    
                    if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                        pp[midx].inb[pp[midx].nnb] = nbas;
                        pp[midx].ipp[pp[midx].nnb] = ncand;
                        pp[midx].zpp[pp[midx].nnb] = M[ncand];
                        (pp[midx].nnb)++; MemCheck(pp,midx);
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            pp[nmidx].inb[pp[nmidx].nnb] = bas;
                            pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                            pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                            (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                        }
                    } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                        pp[midx].inb[pploc] = nbas;
                        pp[midx].ipp[pploc] = ncand;
                        pp[midx].zpp[pploc] = M[ncand];
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            // search for the position of the existing PP in the neighbor basin's list
                            c = 0;
                            while (pp[nmidx].inb[c] != bas) { c++; }
                            pp[nmidx].inb[c] = bas;
                            pp[nmidx].ipp[c] = ncand;
                            pp[nmidx].zpp[c] = M[ncand];
                        }
                    }
                }
            }
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
    
    bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to
    
    if (!sinks[bas]) { /* weights of sinks remain zero */
        
        for (k=0; k<8; k++) {
            if (!skip[k]) {
                p = i + di[k];
                q = j + dj[k];
                
                nbas = (int)Basin[K*q+p];
                
                if (nbas != bas) { // if the neighbor is part of a different basin
                    if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                        ncand = K*q+p; // the neighbor is a candidate pour point
                    } else { // the neighbor has equal or lower elevation
                        ncand = K*j+i; // the present point is a candidate pour point
                    }
                    
                    // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                    midx = (int)minima[bas]; // the list position of the minimum that this basin drains to
                    
                    lowestppyet = 0;
                    ppexists = 0;
                    for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                        if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                            ppexists = 1;
                            pploc = c;
                            if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                                lowestppyet = 1;
                            }
                            break;
                        }
                    }
                    
                    if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                        pp[midx].inb[pp[midx].nnb] = nbas;
                        pp[midx].ipp[pp[midx].nnb] = ncand;
                        pp[midx].zpp[pp[midx].nnb] = M[ncand];
                        (pp[midx].nnb)++; MemCheck(pp,midx);
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            pp[nmidx].inb[pp[nmidx].nnb] = bas;
                            pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                            pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                            (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                        }
                    } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                        pp[midx].inb[pploc] = nbas;
                        pp[midx].ipp[pploc] = ncand;
                        pp[midx].zpp[pploc] = M[ncand];
                        
                        // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                        if (!sinks[nbas]) {
                            nmidx = (int)minima[nbas];
                            // search for the position of the existing PP in the neighbor basin's list
                            c = 0;
                            while (pp[nmidx].inb[c] != bas) { c++; }
                            pp[nmidx].inb[c] = bas;
                            pp[nmidx].ipp[c] = ncand;
                            pp[nmidx].zpp[c] = M[ncand];
                        }
                    }
                }
            }
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

bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to

if (!sinks[bas]) { /* weights of sinks remain zero */

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            p = i + di[k];
            q = j + dj[k];

            nbas = (int)Basin[K*q+p];

            if (nbas != bas) { // if the neighbor is part of a different basin
                if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                    ncand = K*q+p; // the neighbor is a candidate pour point
                } else { // the neighbor has equal or lower elevation
                    ncand = K*j+i; // the present point is a candidate pour point
                }

                // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                midx = (int)minima[bas]; // the list position of the minimum that this basin drains to

                lowestppyet = 0;
                ppexists = 0;
                for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                    if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                        ppexists = 1;
                        pploc = c;
                        if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                            lowestppyet = 1;
                        }
                        break;
                    }
                }

                if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                    pp[midx].inb[pp[midx].nnb] = nbas;
                    pp[midx].ipp[pp[midx].nnb] = ncand;
                    pp[midx].zpp[pp[midx].nnb] = M[ncand];
                    (pp[midx].nnb)++; MemCheck(pp,midx);

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        pp[nmidx].inb[pp[nmidx].nnb] = bas;
                        pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                        pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                        (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                    }
                } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                    pp[midx].inb[pploc] = nbas;
                    pp[midx].ipp[pploc] = ncand;
                    pp[midx].zpp[pploc] = M[ncand];

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        // search for the position of the existing PP in the neighbor basin's list
                        c = 0;
                        while (pp[nmidx].inb[c] != bas) { c++; }
                        pp[nmidx].inb[c] = bas;
                        pp[nmidx].ipp[c] = ncand;
                        pp[nmidx].zpp[c] = M[ncand];
                    }
                }
            }
        }
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

bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to

if (!sinks[bas]) { /* weights of sinks remain zero */

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            p = i + di[k];
            q = j + dj[k];

            nbas = (int)Basin[K*q+p];

            if (nbas != bas) { // if the neighbor is part of a different basin
                if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                    ncand = K*q+p; // the neighbor is a candidate pour point
                } else { // the neighbor has equal or lower elevation
                    ncand = K*j+i; // the present point is a candidate pour point
                }

                // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                midx = (int)minima[bas]; // the list position of the minimum that this basin drains to

                lowestppyet = 0;
                ppexists = 0;
                for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                    if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                        ppexists = 1;
                        pploc = c;
                        if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                            lowestppyet = 1;
                        }
                        break;
                    }
                }

                if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                    pp[midx].inb[pp[midx].nnb] = nbas;
                    pp[midx].ipp[pp[midx].nnb] = ncand;
                    pp[midx].zpp[pp[midx].nnb] = M[ncand];
                    (pp[midx].nnb)++; MemCheck(pp,midx);

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        pp[nmidx].inb[pp[nmidx].nnb] = bas;
                        pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                        pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                        (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                    }
                } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                    pp[midx].inb[pploc] = nbas;
                    pp[midx].ipp[pploc] = ncand;
                    pp[midx].zpp[pploc] = M[ncand];

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        // search for the position of the existing PP in the neighbor basin's list
                        c = 0;
                        while (pp[nmidx].inb[c] != bas) { c++; }
                        pp[nmidx].inb[c] = bas;
                        pp[nmidx].ipp[c] = ncand;
                        pp[nmidx].zpp[c] = M[ncand];
                    }
                }
            }
        }
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

bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to

if (!sinks[bas]) { /* weights of sinks remain zero */

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            p = i + di[k];
            q = j + dj[k];

            nbas = (int)Basin[K*q+p];

            if (nbas != bas) { // if the neighbor is part of a different basin
                if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                    ncand = K*q+p; // the neighbor is a candidate pour point
                } else { // the neighbor has equal or lower elevation
                    ncand = K*j+i; // the present point is a candidate pour point
                }

                // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                midx = (int)minima[bas]; // the list position of the minimum that this basin drains to

                lowestppyet = 0;
                ppexists = 0;
                for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                    if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                        ppexists = 1;
                        pploc = c;
                        if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                            lowestppyet = 1;
                        }
                        break;
                    }
                }

                if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                    pp[midx].inb[pp[midx].nnb] = nbas;
                    pp[midx].ipp[pp[midx].nnb] = ncand;
                    pp[midx].zpp[pp[midx].nnb] = M[ncand];
                    (pp[midx].nnb)++; MemCheck(pp,midx);

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        pp[nmidx].inb[pp[nmidx].nnb] = bas;
                        pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                        pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                        (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                    }
                } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                    pp[midx].inb[pploc] = nbas;
                    pp[midx].ipp[pploc] = ncand;
                    pp[midx].zpp[pploc] = M[ncand];

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        // search for the position of the existing PP in the neighbor basin's list
                        c = 0;
                        while (pp[nmidx].inb[c] != bas) { c++; }
                        pp[nmidx].inb[c] = bas;
                        pp[nmidx].ipp[c] = ncand;
                        pp[nmidx].zpp[c] = M[ncand];
                    }
                }
            }
        }
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

bas = (int)Basin[K*j+i]; // the basin that this point drains to, which is also the index of the point that it ultimately drains to

if (!sinks[bas]) { /* weights of sinks remain zero */

    for (k=0; k<8; k++) {
        if (!skip[k]) {
            p = i + di[k];
            q = j + dj[k];

            nbas = (int)Basin[K*q+p];

            if (nbas != bas) { // if the neighbor is part of a different basin
                if (M[K*q+p] > M[K*j+i]) { // if the neighbor has higher elevation
                    ncand = K*q+p; // the neighbor is a candidate pour point
                } else { // the neighbor has equal or lower elevation
                    ncand = K*j+i; // the present point is a candidate pour point
                }

                // compare the candidate point to existing PP to this neighbor basin, if it exists. If it's lower or there isn't one yet, make it the PP to this neighbor
                midx = (int)minima[bas]; // the list position of the minimum that this basin drains to

                lowestppyet = 0;
                ppexists = 0;
                for (c=0; c<pp[midx].nnb; c++) { // the loop shouldn't even execute once if there are no pp's yet for this basin (pp[midx].nnb == 0)
                    if (pp[midx].inb[c] == nbas) { // if this is an existing entry for the neighboring basin in question
                        ppexists = 1;
                        pploc = c;
                        if (M[ncand] < pp[midx].zpp[c]) { // if the candidate PP is lower than the existing one
                            lowestppyet = 1;
                        }
                        break;
                    }
                }

                if (!ppexists) { // if no pp exists yet for this pair of basins, record the candidate as the pp
                    pp[midx].inb[pp[midx].nnb] = nbas;
                    pp[midx].ipp[pp[midx].nnb] = ncand;
                    pp[midx].zpp[pp[midx].nnb] = M[ncand];
                    (pp[midx].nnb)++; MemCheck(pp,midx);

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        pp[nmidx].inb[pp[nmidx].nnb] = bas;
                        pp[nmidx].ipp[pp[nmidx].nnb] = ncand;
                        pp[nmidx].zpp[pp[nmidx].nnb] = M[ncand];
                        (pp[nmidx].nnb)++; MemCheck(pp,nmidx);
                    }
                } else if (lowestppyet) { // if there was an existing PP, but the candidate is lower, replace the existing one
                    pp[midx].inb[pploc] = nbas;
                    pp[midx].ipp[pploc] = ncand;
                    pp[midx].zpp[pploc] = M[ncand];

                    // also record this as the PP from the neighbor basin to this one, provided the neighbor isn't a sink basin
                    if (!sinks[nbas]) {
                        nmidx = (int)minima[nbas];
                        // search for the position of the existing PP in the neighbor basin's list
                        c = 0;
                        while (pp[nmidx].inb[c] != bas) { c++; }
                        pp[nmidx].inb[c] = bas;
                        pp[nmidx].ipp[c] = ncand;
                        pp[nmidx].zpp[c] = M[ncand];
                    }
                }
            }
        }
    }
}
    

    
} // end GetPourPoints



void TracePaths(pourpoint pp[],const int K, const int J, double M[], double Basin[], double minima[], double minimaIdx[], int nMin, double s[], double b[])
{

int i, j, idx, n, c, cmin, nc, g, wing, nbas, winner, loser, reachedsink, loopedback, thisbaslist, nextbaslist, nextbasidx, ppexists, pploc;
double zppmin;    

// printf("nMin=%d\n",nMin);

// For each basin, trace the path of lowest PPs.
for (n=0; n<nMin; n++) {
    
    if (minimaIdx[n] != -1) { // check to make sure this minimum hasn't been merged into another basin (if it has, it doesn't exist anymore).
        reachedsink = 0;
        loopedback = 0;
        nextbasidx = (int)minimaIdx[n];
        
//         printf("Tracing n=%d idx=%d: ",n,(int)minimaIdx[n]);
        
        while (!reachedsink && !loopedback) {
            
            thisbaslist = (int)minima[nextbasidx];
            nextbasidx = pp[thisbaslist].inb[pp[thisbaslist].lpp];
            nextbaslist = (int)minima[nextbasidx];
            
//             printf("->%d",(int)minimaIdx[thisbaslist]);
            
            if (s[nextbasidx]) { // if the next basin is a sink
//                 printf("->SINK %d\n",nextbasidx);
                reachedsink = 1;
            } else if (pp[nextbaslist].inb[pp[nextbaslist].lpp] == (int)minimaIdx[thisbaslist]) { // if the next basin drains back to this one
//                 printf("->LOOP %d->%d\n",nextbasidx,pp[nextbaslist].inb[pp[nextbaslist].lpp]);
                loopedback = 1;
            } // otherwise we've not reached a sink or looped back; keep tracing the path
        }
        if (loopedback) { // If we looped back to the starting basin, aggregate basins, find lowest PP, and trace lowest PP path again
            
            // of the two basins, find the one with the lowest PP that does not connect it to the other. This is the new identity for the aggregate basin
            zppmin = INF;
            for (c=0; c<pp[thisbaslist].nnb; c++) {
                if (c != pp[thisbaslist].lpp && pp[thisbaslist].zpp[c] < zppmin) { // look for lowest pour point, excluding the original
                    cmin = c;
                    winner = thisbaslist;
                    loser = nextbaslist;
                    zppmin = pp[winner].zpp[c];
                }
            }
            for (c=0; c<pp[nextbaslist].nnb; c++) {
                if (c != pp[nextbaslist].lpp && pp[nextbaslist].zpp[c] < zppmin) { // look for lowest pour point, excluding the original
                    cmin = c;
                    winner = nextbaslist;
                    loser = thisbaslist;
                    zppmin = pp[winner].zpp[c];
                }
            }
            
            // Delete the entry for the loser from the winner's table
            for (c=pp[winner].lpp; c<pp[winner].nnb; c++) {
                pp[winner].inb[c] = pp[winner].inb[c+1];
                pp[winner].ipp[c] = pp[winner].ipp[c+1];
                pp[winner].zpp[c] = pp[winner].zpp[c+1];
            }
            (pp[winner].nnb)--;
            
            // mark the new lowest PP for the winner
            if (cmin > pp[winner].lpp) { // if the new lowest PP was further down the list than the loser (which is now gone from the list), the relevant entry is now one higher on the list
                pp[winner].lpp = cmin-1;
            } else {
                pp[winner].lpp = cmin;
            }
            
            // visit each of the loser's "partner" basins with which it shares a PP, and deal with those PPs, now that the loser and winner are a single basin
            for (c=0; c<pp[loser].nnb; c++) {
                if (c != pp[loser].lpp) { // if it's not the original pp (which doesn't point back here anymore)
                    
                    nbas = (int)minima[pp[loser].inb[c]]; // list position of the neighbor basin
                    
                    // find whether the winner already shares a PP with this basin, and if so, record the position in the winner's list
                    ppexists = 0;
                    nc = 0;
                    while (!ppexists && nc<pp[winner].nnb) {
                        if (pp[winner].inb[nc] == pp[loser].inb[c]) {
                            ppexists = 1;
                            pploc = nc;
                        }
                        nc++;
                    }
                    
                    
                    if (!s[pp[loser].inb[c]]) { // if the neighbor is not a sink (which don't have recorded PPs)
                        
                        nc = 0;
                        while (pp[nbas].inb[nc] != (int)minimaIdx[loser]) { nc++; } // Find the list position of the neighbor's PP that points to the loser
                        
                        if (ppexists) { // if there was already a shared PP
                            if (pp[loser].zpp[c] < pp[winner].zpp[pploc]) { // if the one from the loser to the loser-neighbor is lower than the one from the winner to the loser-neighbor
                                // replace the existing entry in the winner's list
                                pp[winner].inb[pploc] = pp[loser].inb[c];
                                pp[winner].ipp[pploc] = pp[loser].ipp[c];
                                pp[winner].zpp[pploc] = pp[loser].zpp[c];
                                
                                // find the existing entry in the loser-neighbor's list that pointed to the winner
                                wing = 0;
                                while (pp[nbas].inb[wing] != (int)minimaIdx[winner]) { wing++; }
                                
                                // make the entry in the loser-neighbor's list that pointed to the loser now point to the winner instead
                                pp[nbas].inb[nc] = (int)minimaIdx[winner];
                                
                                // delete the existing entry in the loser-neighbor's list that pointed to the winner
                                for (g=wing; g<pp[nbas].nnb; g++) {
                                    pp[nbas].inb[g] = pp[nbas].inb[g+1];
                                    pp[nbas].ipp[g] = pp[nbas].ipp[g+1];
                                    pp[nbas].zpp[g] = pp[nbas].zpp[g+1];
                                }
                                (pp[nbas].nnb)--;
                                if (pp[nbas].lpp > wing) { // if the loser-neighbor's lowest PP was further down the list than the winner
                                    (pp[nbas].lpp)--;
                                } else if (pp[nbas].lpp == wing) { // if the loser-neighbor's lowest PP WAS the winner (this shouldn't happen, because we arrived here having found that the PP to the loser was lower)
                                    printf("Please notify the developer of a type 1 error\n");
                                }
                            } else { // if the existing PP is lower than or equal to the PP from the loser-neighbor to the loser
                                // do nothing to the winner's list
                                
                                // delete the entry in the loser-neighbor's list that points to the loser
                                for (g=nc; g<pp[nbas].nnb; g++) {
                                    pp[nbas].inb[g] = pp[nbas].inb[g+1];
                                    pp[nbas].ipp[g] = pp[nbas].ipp[g+1];
                                    pp[nbas].zpp[g] = pp[nbas].zpp[g+1];
                                }
                                (pp[nbas].nnb)--;
                                if (pp[nbas].lpp > nc) { // if the loser-neighbor's lowest PP was further down the list than the loser (which is now gone from the list), the relevant entry is now one higher on the list
                                    (pp[nbas].lpp)--;
                                } else if (pp[nbas].lpp == nc) { // If this is true (loser-neighbor's lowest PP WAS the loser), it must have been tied with the PP from the loser-neighbor to the winner.
                                    // find the entry in the loser-neighbor's list that points to the winner, and make it the new lowest PP
                                    wing = 0;
                                    while (pp[nbas].inb[wing] != (int)minimaIdx[winner]) { wing++; }
                                    pp[nbas].lpp = wing;
                                }
                            }
                        } else { // if there was no existing shared PP
                            // add the entry to the winner's list pointing to the loser-neighbor
                            pp[winner].inb[pp[winner].nnb] = pp[loser].inb[c];
                            pp[winner].ipp[pp[winner].nnb] = pp[loser].ipp[c];
                            pp[winner].zpp[pp[winner].nnb] = pp[loser].zpp[c];
                            (pp[winner].nnb)++; // increment the number of neighbor basins in the list, since we just added one
                            // note that there is no possibility that we just added the winner's lowest PP, because by definition, the winner's lowest PP was lower than all the loser's PPs
                            MemCheck(pp,winner); // Make sure we haven't exceeded the number of allowed neighboring basins

                            
                            // make the entry in the loser-neighbor's list point to the winner instead of the loser
                            pp[nbas].inb[nc] = (int)minimaIdx[winner];
                        }
                    } else { // the neighbor is a sink
                        if (ppexists) { // if there was already a shared PP
                            if (pp[loser].zpp[c] < pp[winner].zpp[pploc]) { // if the one from the loser to the loser-neighbor is lower
                                // replace the old with the new in the winner's PP list
                                pp[winner].inb[pploc] = pp[loser].inb[c];
                                pp[winner].ipp[pploc] = pp[loser].ipp[c];
                                pp[winner].zpp[pploc] = pp[loser].zpp[c];
                            } // otherwise  do nothing, because we need not modify the winner's list, and there is no loser-neighbor's list
                            // Note that there is no possibility that we just replaced the winner's lowest PP, because by definition, the winner's lowest PP was lower than all the loser's PPs
                        } else { // if there was no existing shared PP
                            // add the entry to the winner's list pointing to the loser-neighbor-sink
                            pp[winner].inb[pp[winner].nnb] = pp[loser].inb[c];
                            pp[winner].ipp[pp[winner].nnb] = pp[loser].ipp[c];
                            pp[winner].zpp[pp[winner].nnb] = pp[loser].zpp[c];
                            (pp[winner].nnb)++; // increment the number of neighbor basins in the list, since we just added one
                            MemCheck(pp,winner); // Make sure we haven't exceeded the number of allowed neighboring basins

                        }
                    }
                }
            }
            
            
            // "repaint" the loser with the index of the minimum of the winner by calling "Repaint", a function similar to GetBasin that just expands out recursively by basin number.
            idx = (int)minimaIdx[loser];
            i = idx % K;
            j = (idx - i)/K;
            
            Repaint(i, j, idx, (int)minimaIdx[winner], K, J, Basin, b);
            
            // minima[(int)minimaIdx[loser]] = 0; // The loser basin no longer exists.
            minimaIdx[loser] = -1; // This allows us to check at the beginning of this loop to make sure basin n is still a valid basin
               
            // we didn't make it to a sink; decrement the counter so once we exit this iteration of the loop we'll try again. If n was the winner, we'll see if we make it to a sink. If n was the loser, we'll skip it, and come back to the winner, which has higher n, in a later iteration.
            // A note on efficiency: really we only need to start re-tracing from the winner, because everything upstream is unchanged. As long as we keep track of where we are in the nMin queue...
            n--;
        } // else if (reachedsink) { } // If we reached a sink basin (In this version, we don't do anything here)
    }
        
} // finished tracing paths

} // end TracePaths



void DrainToMe(const int i, const int j, const int K, const int J, double label, double Z, double M[], double D[], double F[], double Basin[], double bdy[])
{

int k, p, q, bl, br, bu, bd;
int r[] = {5,6,7,8,1,2,3,4};
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};


Basin[K*j+i] = -1; // this point has now been visited


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


if (i==0) {
    
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

} else if (i==K-1) {

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


}

if (j==0) {
    
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
    
} else if (j==J-1) {

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

}


// loop through neighbors. To make the drainage directions in the basin pretty, we'll have to come up with a way of visiting neighbors in order of distance from minimum or something.
for (k=0; k<8; k++) {

    if (!skip[k]) {
        
        p = i + di[k];
        q = j + dj[k];
        
        if (Basin[K*q+p] == label && M[K*q+p] <= Z) { // if neighbor k is in the basin, has not yet been visited, and is at or below the flooding level
            D[K*q+p] = r[k]; // make neighbor k drain to i,j.
            M[K*q+p] = Z; // raise elevation to the water level
            F[K*q+p] = 1; // it's been flooded
            DrainToMe(p,q,K,J,label,Z,M,D,F,Basin,bdy); // call DrainToMe recursively
        }        
    }
}

} // end DrainToMe()


void Flood(const int K, const int J, double M[], double Mfl[], double D[], double F[], double Basin[], double minima[], double minimaIdx[], int *numMin, double *a, double s[], double b[]) 
{
    
    int i, j, n, c, idx, nMin;
    double zppmin;
    pourpoint *pp;

    D8Dir(K, J, M, D, minima, minimaIdx, numMin, a, b, s); // Calculate D8 drainage directions, recording minima
    
    
    /* find the D8 basin that drains to the minimum. */
    nMin = *numMin; // Once everything is working, it won't be necessary to decrement numMin to make sure we're taking care of each minimum, so we won't need to pass it to functions below
//     pourpoint pp[nMin];
    pp = mxCalloc(nMin, sizeof(pourpoint)); // allocate memory for pour point table

    for (n=0; n<nMin; n++) {
        
        pp[n].nnb = 0; // initialize number of neighboring basins for each min to zero
        
        // find next minimum
        idx = (int)minimaIdx[n];
        i = idx % K;
        j = (idx - i)/K;
        GetBasin(i, j, idx, K, J, D, Basin, b);
    }
    
    // Now mark basins that drain to sinks
    for (i=0; i<K; i++) {
        for (j=0; j<J; j++) {
            if (s[K*j+i]) {
                GetBasin(i, j, K*j+i, K, J, D, Basin, b);
            }
        }
    }
    
    
    // NOW FIND POUR POINTS, passing an array of structs
    GetPourPoints(pp, K, J, M, Basin, minima, minimaIdx, numMin, s, b);
    
    // NOW LOOP THROUGH MINIMA AGAIN AND FIND THE LOWEST POUR POINT FOR EACH BASIN
    for (n=0; n<nMin; n++) {
        zppmin = INF;
        for (c=0; c<pp[n].nnb; c++) {
            if (pp[n].zpp[c] < zppmin) {
                pp[n].lpp = c;
                zppmin = pp[n].zpp[c];
            }
        }
    }
    
    TracePaths(pp, K, J, M, Basin, minima, minimaIdx, nMin, s, b);
    
    // loop through minima, flooding each one to the level of the lowest PP
    for (n=0; n<nMin; n++) {
        if (minimaIdx[n] != -1) {
            idx = pp[n].ipp[pp[n].lpp];
            i = idx % K;
            j = (idx - i)/K;
            DrainToMe(i, j, K, J, minimaIdx[n], pp[n].zpp[pp[n].lpp], Mfl, D, F, Basin, b);
        }
    }
    
    // deallocate memory
    mxFree(pp);
    
} // end Flood



void GetWeights(const int K, const int J, double M[], double W[], double minima[], double *a, int *numMin, double bdy[], double sinks[])
{

int i, j, k, Kj;
int bl, br, bu, bd;
int downslope[8];
int skip[] = {0,0,0,0,0,0,0,0};
int di[8], di0[] = {0,-1,-1,-1,0,1,1,1};
int dj[8], dj0[] = {1,1,0,-1,-1,-1,0,1};
double b[8], s, sumS, invsumS, e0, en;


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


// set b, the inverse distance to each neighbor, with distance normalized to 1/dx
for (k=0; k<8; k++) {
    switch (k) {
        case 0:
        case 4:
            b[k] = 1;
            break;
        case 2:
        case 6:
            b[k] = 1/(*a);
            break;
        case 1:
        case 3:
        case 5:
        case 7:
            b[k] = 1/sqrt(1 + (*a)*(*a));
            break;
    }
}


/* INTERIOR POINTS */

for (k=0; k<8; k++) {di[k]=di0[k];}
for (k=0; k<8; k++) {dj[k]=dj0[k];}

for (j=1; j<(J-1); j++) { // Loop through columns
    Kj=K*j;
    for (i=1; i<(K-1); i++) { // Loop through rows
        sumS = 0; // Set sum of slopes to zero        
        if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink            
            e0=M[Kj+i]; // Elevation at (i,j)            
            for (k=0; k<8; k++) { // Loop through neighbors                
                downslope[k] = 0; // reset the downslope indicator                
                en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k                
                // Only proceed to investigate this facet if it is a downslope facet
                if (en < e0) {                    
                    downslope[k] = 1;
                    s = (e0-en)*b[k]; // slope to neighbor k
                    W[k*J*K+Kj+i] = s; 
                    sumS += s;                    
                }
            }
        }
        
        // If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
        if (sumS) {            
            invsumS = 1/sumS;            
            for (k=0; k<8; k++) { // Loop through facets
                if (downslope[k]) {                    
                    W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1                    
                }
            }            
        } else { // It is a local minimum (this includes sinks)
            minima[Kj+i] = 1;
            (*numMin)++;
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
    skip[2]=1;
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

for (i=1; i<(K-1); i++) { // Loop through rows
    sumS = 0; // Set sum of slopes to zero
    if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
        e0=M[Kj+i]; // Elevation at (i,j)
        for (k=0; k<8; k++) { // Loop through neighbors
            downslope[k] = 0; // reset the downslope indicator
            if (!skip[k]) {
                en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
                // Only proceed to investigate this facet if it is a downslope facet
                if (en < e0) {
                    downslope[k] = 1;
                    s = (e0-en)*b[k]; // slope to neighbor k
                    W[k*J*K+Kj+i] = s;
                    sumS += s;
                }
            }
        }
    }
    
    // If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
    if (sumS) {
        invsumS = 1/sumS;
        for (k=0; k<8; k++) { // Loop through facets
            if (downslope[k]) {
                W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
            }
        }
    } else { // It is a local minimum (this includes sinks)
        minima[Kj+i] = 1;
        (*numMin)++;
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
    skip[6]=1;
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

for (i=1; i<(K-1); i++) { // Loop through rows
    sumS = 0; // Set sum of slopes to zero
    if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
        e0=M[Kj+i]; // Elevation at (i,j)
        for (k=0; k<8; k++) { // Loop through neighbors
            downslope[k] = 0; // reset the downslope indicator
            if (!skip[k]) {
                en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
                // Only proceed to investigate this facet if it is a downslope facet
                if (en < e0) {
                    downslope[k] = 1;
                    s = (e0-en)*b[k]; // slope to neighbor k
                    W[k*J*K+Kj+i] = s;
                    sumS += s;
                }
            }
        }
    }
    
    // If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
    if (sumS) {
        invsumS = 1/sumS;
        for (k=0; k<8; k++) { // Loop through facets
            if (downslope[k]) {
                W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
            }
        }
    } else { // It is a local minimum (this includes sinks)
        minima[Kj+i] = 1;
        (*numMin)++;
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
    skip[0]=1;
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

for (j=1; j<(J-1); j++) { // Loop through rows
    Kj = K*j;
    sumS = 0; // Set sum of slopes to zero
    if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
        e0=M[Kj+i]; // Elevation at (i,j)
        for (k=0; k<8; k++) { // Loop through neighbors
            downslope[k] = 0; // reset the downslope indicator
            if (!skip[k]) {
                en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
                // Only proceed to investigate this facet if it is a downslope facet
                if (en < e0) {
                    downslope[k] = 1;
                    s = (e0-en)*b[k]; // slope to neighbor k
                    W[k*J*K+Kj+i] = s;
                    sumS += s;
                }
            }
        }
    }
    
    // If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
    if (sumS) {
        invsumS = 1/sumS;
        for (k=0; k<8; k++) { // Loop through facets
            if (downslope[k]) {
                W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
            }
        }
    } else { // It is a local minimum (this includes sinks)
        minima[Kj+i] = 1;
        (*numMin)++;
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
    skip[4]=1;
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

for (j=1; j<(J-1); j++) { // Loop through rows
    Kj = K*j;
    sumS = 0; // Set sum of slopes to zero
    if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
        e0=M[Kj+i]; // Elevation at (i,j)
        for (k=0; k<8; k++) { // Loop through neighbors
            downslope[k] = 0; // reset the downslope indicator
            if (!skip[k]) {
                en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
                // Only proceed to investigate this facet if it is a downslope facet
                if (en < e0) {
                    downslope[k] = 1;
                    s = (e0-en)*b[k]; // slope to neighbor k
                    W[k*J*K+Kj+i] = s;
                    sumS += s;
                }
            }
        }
    }
    
    // If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
    if (sumS) {
        invsumS = 1/sumS;
        for (k=0; k<8; k++) { // Loop through facets
            if (downslope[k]) {
                W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
            }
        }
    } else { // It is a local minimum (this includes sinks)
        minima[Kj+i] = 1;
        (*numMin)++;
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
    skip[2]=1;
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
    skip[0]=1;
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

sumS = 0; // Set sum of slopes to zero
if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
    e0=M[Kj+i]; // Elevation at (i,j)
    for (k=0; k<8; k++) { // Loop through neighbors
        downslope[k] = 0; // reset the downslope indicator
        if (!skip[k]) {
            en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
            // Only proceed to investigate this facet if it is a downslope facet
            if (en < e0) {
                downslope[k] = 1;
                s = (e0-en)*b[k]; // slope to neighbor k
                W[k*J*K+Kj+i] = s;
                sumS += s;
            }
        }
    }
}

// If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
if (sumS) {
    invsumS = 1/sumS;
    for (k=0; k<8; k++) { // Loop through facets
        if (downslope[k]) {
            W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
        }
    }
} else { // It is a local minimum (this includes sinks)
    minima[Kj+i] = 1;
    (*numMin)++;
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
    skip[6]=1;
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
    skip[0]=1;
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

sumS = 0; // Set sum of slopes to zero
if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
    e0=M[Kj+i]; // Elevation at (i,j)
    for (k=0; k<8; k++) { // Loop through neighbors
        downslope[k] = 0; // reset the downslope indicator
        if (!skip[k]) {
            en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
            // Only proceed to investigate this facet if it is a downslope facet
            if (en < e0) {
                downslope[k] = 1;
                s = (e0-en)*b[k]; // slope to neighbor k
                W[k*J*K+Kj+i] = s;
                sumS += s;
            }
        }
    }
}

// If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
if (sumS) {
    invsumS = 1/sumS;
    for (k=0; k<8; k++) { // Loop through facets
        if (downslope[k]) {
            W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
        }
    }
} else { // It is a local minimum (this includes sinks)
    minima[Kj+i] = 1;
    (*numMin)++;
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
    skip[2]=1;
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
    skip[4]=1;
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

sumS = 0; // Set sum of slopes to zero
if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
    e0=M[Kj+i]; // Elevation at (i,j)
    for (k=0; k<8; k++) { // Loop through neighbors
        downslope[k] = 0; // reset the downslope indicator
        if (!skip[k]) {
            en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
            // Only proceed to investigate this facet if it is a downslope facet
            if (en < e0) {
                downslope[k] = 1;
                s = (e0-en)*b[k]; // slope to neighbor k
                W[k*J*K+Kj+i] = s;
                sumS += s;
            }
        }
    }
}

// If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
if (sumS) {
    invsumS = 1/sumS;
    for (k=0; k<8; k++) { // Loop through facets
        if (downslope[k]) {
            W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
        }
    }
} else { // It is a local minimum (this includes sinks)
    minima[Kj+i] = 1;
    (*numMin)++;
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
    skip[6]=1;
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
    skip[4]=1;
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

sumS = 0; // Set sum of slopes to zero
if (!sinks[Kj+i]) { // Only calculate weights if (i,j) is not a user-specified sink
    e0=M[Kj+i]; // Elevation at (i,j)
    for (k=0; k<8; k++) { // Loop through neighbors
        downslope[k] = 0; // reset the downslope indicator
        if (!skip[k]) {
            en = M[K*(j+dj[k])+(i+di[k])]; // elevation of neighbor k
            // Only proceed to investigate this facet if it is a downslope facet
            if (en < e0) {
                downslope[k] = 1;
                s = (e0-en)*b[k]; // slope to neighbor k
                W[k*J*K+Kj+i] = s;
                sumS += s;
            }
        }
    }
}

// If (i,j) is not a local minimum or a sink (sumS remains zero), assign weights; otherwise, W[i,j,:] remains zero and the location is flagged in minima[]
if (sumS) {
    invsumS = 1/sumS;
    for (k=0; k<8; k++) { // Loop through facets
        if (downslope[k]) {
            W[k*J*K+Kj+i] = invsumS*W[k*J*K+Kj+i]; // make weights sum to 1
        }
    }
} else { // It is a local minimum (this includes sinks)
    minima[Kj+i] = 1;
    (*numMin)++;
}



} /* end GetWeights() */




void FixMirrorWeights(const int K, const int J, double W[], double bdy[]) {
    
int i, j, k1, k2, k3;

/* left */
if (bdy[0] == 1) {
    j = 0;
    k1 = 7;
    k2 = 0;
    k3 = 1;
    
    for (i=1; i<(K-1); i++) {
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
    
    for (i=1; i<(K-1); i++) {
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
    
    for (j=1; j<(J-1); j++) {
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
    
    for (j=1; j<(J-1); j++) {
        W[k1*K*J+K*j+i] *= 0.5;
        W[k2*K*J+K*j+i] *= 0.5;
        W[k3*K*J+K*j+i] *= 0.5;
    }
}
    
} /* end FixMirrorWeights() */


void ReplaceWeights(const int K, const int J, double W[], double D[], double F[]) {

    int i, j, k;

    for (i=0; i<K; i++) {
        for (j=0; j<J; j++) {
            if (F[K*j+i]) {
                for (k=0; k<8; k++) {
                    W[k*K*J+K*j+i]=0; // Wipe it clean
                }
                W[(int) (D[K*j+i]-1)*K*J+K*j+i]=1; // Replace it with the D8 weight
            }
        }
    }

} // end ReplaceWeights()



void GetArea(const int i, const int j, const int K, const int J, double A[], double W[], double bdy[])
{

int k, p, q, Kqplusp, c=K*j+i, widx, bl, br, bu, bd;
int r[] = {4,5,6,7,0,1,2,3};
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};
double theW, w[] = {1,1,1,1,1,1,1,1};

if (A[c]>0) {return;} // Bail out if we already know the area

A[c] = 1; /* each cell drains at least itself */

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


if (i==0) {
    
    switch (bu) {
    case FIXED: // fixed upper
        skip[1]=1;
        skip[2]=1;
        skip[3]=1;
        break;
    case MIRROR: // mirror upper
        skip[1]=1;
        skip[2]=1;
        skip[3]=1;
        w[5] = 2;
        w[6] = 2;
        w[7] = 2;
        break;
    case PERIODIC: // periodic upper
        di[1]=K-1;
        di[2]=K-1;
        di[3]=K-1;
        break;
    }

} else if (i==K-1) {

    switch (bd) {
    case FIXED: // fixed lower
        skip[5]=1;
        skip[6]=1;
        skip[7]=1;
        break;
    case MIRROR: // mirror lower
        skip[5]=1;
        skip[6]=1;
        skip[7]=1;
        w[1] = 2;
        w[2] = 2;
        w[3] = 2;
        break;
    case PERIODIC: // periodic lower
        di[5]=1-K;
        di[6]=1-K;
        di[7]=1-K;
        break;
    }


}

if (j==0) {
    
    switch (bl) {
    case FIXED: // fixed left
        skip[3]=1;
        skip[4]=1;
        skip[5]=1;
        break;
    case MIRROR: // mirror left
        skip[3]=1;
        skip[4]=1;
        skip[5]=1;
        w[7] = 2;
        w[0] = 2;
        w[1] = 2;
        break;
    case PERIODIC: // periodic left
        dj[3]=J-1;
        dj[4]=J-1;
        dj[5]=J-1;
        break;
    }
    
} else if (j==J-1) {

    switch (br) {
    case FIXED: // fixed right
        skip[7]=1;
        skip[0]=1;
        skip[1]=1;
        break;
    case MIRROR: // mirror right
        skip[7]=1;
        skip[0]=1;
        skip[1]=1;
        w[3] = 2;
        w[4] = 2;
        w[5] = 2;
        break;
    case PERIODIC: // periodic right
        dj[7]=1-J;
        dj[0]=1-J;
        dj[1]=1-J;
        break;
    }

}


for (k=0; k<8; k++) {  /* loop through drainage directions */
    
    if (!skip[k]) {
        
        p = i + di[k];
        q = j + dj[k];
        
        Kqplusp = K*q+p;
        widx=r[k]*J*K+Kqplusp;
        
        if ((theW=W[widx])>0) {  /* if the current drainage direction has a weight */
            
            GetArea(p, q, K, J, A, W, bdy);  /* recursive call to get drainage area for neighbor k */
            A[c] += w[k] * theW*A[Kqplusp]; /* increment A(i,j) by the weighted area of k */
            /* Purpose of the w prefactor: A point on a mirrored boundary receiving flow from a point off the boundary should receive twice the area */   
        }
    }
    
} /* end of drainage direction loop */

} /* end GetArea subfunction */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    mxArray *Wptr, *Mptr, *Dptr, *Basinptr, *minimaIdxptr;
    double *M, *Mfl, *Basin, *minimaIdx, *D, *F, *W, *A, *minima, *a, *fl, *b, *s, zppmin;
    int i, j, n, c, idx, K, J, nMin, *numMin, ndims=3, dims[]={0, 0, 8};
    
    // Get pointers to inputs
    M = (double *)mxGetPr(prhs[0]); // elevations
    a = (double *)mxGetPr(prhs[1]); // dy/dx
    s = (double *)mxGetPr(prhs[2]); /* sinks where recursive upslope area function should be called */
    b = (double *)mxGetPr(prhs[3]); /* [bleft bright bupper blower] */
    fl = (double *)mxGetPr(prhs[4]); // do or don't do flooding
    
    // Get dimensions of input matrix of elevations
    K=dims[0]=mxGetM(prhs[0]);
    J=dims[1]=mxGetN(prhs[0]);
    
    
    // Create arrays for return arguments
    minima = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // 1 if an element is a local minimum, zero otherwise
    F = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix that indicates what was flooded
    
    numMin = mxCalloc(1, sizeof(int)); // allocate memory for this 1-element integer array
    *numMin = 0;
    
    
    /* FLOODING */
    if (*fl) {  // If user requested flooding
        
        // Create arrays
        D = (double *)mxGetPr(Dptr= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of D8 drainage directions
        minimaIdx = (double *)mxGetPr(minimaIdxptr= mxCreateDoubleMatrix(K, J, mxREAL)); // list of indices of minima
        Basin = (double *)mxGetPr(Basinptr= mxCreateDoubleMatrix(K, J, mxREAL)); /* indicates which cells drain to the selected point (1) and which don't (0) */
        Mfl = (double *)mxGetPr(Mptr= mxCreateDoubleMatrix(K, J, mxREAL)); // the post-landslide elevations

        // copy the input elevations. We do this to avoid modifying the input in memory, which might surprise the user
        for (n=0; n<(K*J); n++) {
            Mfl[n] = M[n];
        }

        // Find the paths that flow would take if closed depressions in the grid overflow
        Flood(K,J,M,Mfl,D,F,Basin,minima,minimaIdx,numMin,a,s,b);
        
        // The weights matrix
        W = (double *)mxGetPr(plhs[3]= mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL));
        
        // get D-infinity weights, using the regraded surface
        GetWeights(K, J, Mfl, W, minima, a, numMin, b, s);
        
        // Replace weights for flooded cells with weights of 1 in the directions indicated by D8
        ReplaceWeights(K, J, W, D, F);
        
        
        /* free memory */
        mxDestroyArray(Dptr);
        mxDestroyArray(Basinptr);
        mxDestroyArray(minimaIdxptr);
        mxDestroyArray(Mptr);        
        
    } else { // flooding was not requested
        
        W = (double *)mxGetPr(plhs[3]= mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL));
        
        // get D-infinity weights
        GetWeights(K, J, M, W, minima, a, numMin, b, s);
        
    }
    
    
    // Free memory
    mxFree(numMin);
    
    /* If necessary, correct weights for points on mirror boundaries */
    /* (A point on a mirrored boundary donating flow to a point off the boundary should only donate half its area) */
    if (b[0]==1 || b[1]==1 || b[2]==1 || b[3]==1) {
        FixMirrorWeights(K, J, W, b);
    }
    
    
    // Create array for return argument A
    A = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); // matrix of TCAs
    
    // Get contributing area for each element
    for (i=0; i<K; i++) { // all points
        for (j=0; j<J; j++) {
            GetArea(i, j, K, J, A, W, b);
        }
    }
    
} // End mexFunction
