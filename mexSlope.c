 /*
 * =============================================================
 * mexSlope.c - C MEX file to compute topographic gradient
 *
 * usage: S = mexSlope(M,dx,dy,bvec) 
          where M is a K x J matrix of
 *        elevations and dx, dy are the cell dimensions.
 *        Returns S, a K x J matrix of topographic slopes computed
 *        as the geometric mean of slope in the 4 directions
 *        defined by a cell and its 8 neighbors.
 *        bvec is a vector of codes for the [left, right, upper, lower]
 *        boundary conditions: 0=fixed, 1=mirror, 2=periodic
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2006 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include <math.h>

#define FIXED    0
#define MIRROR   1
#define PERIODIC 2

/* If necessary, could use this function to pre-screen for local minima */

void GetSlopes(const int K, const int J, double M[], double S[], double *dx, double *dy, double *b)
{

int i, j, Kj, iup, idown, jup, jdown;
int Kjup, Kjdown;
int bl, br, bu, bd;
double invtwodx, invtwody, gx, gy;

// left
if (b[0] == 0) {
    bl=FIXED;
} else if (b[0] == 1) {
    bl=MIRROR;
} else if (b[0] == 2) {
    bl=PERIODIC;
}

// right
if (b[1] == 0) {
    br=FIXED;
} else if (b[1] == 1) {
    br=MIRROR;
} else if (b[1] == 2) {
    br=PERIODIC;
}

// upper
if (b[2] == 0) {
    bu=FIXED;
} else if (b[2] == 1) {
    bu=MIRROR;
} else if (b[2] == 2) {
    bu=PERIODIC;
}

// lower
if (b[3] == 0) {
    bd=FIXED;
} else if (b[3] == 1) {
    bd=MIRROR;
} else if (b[3] == 2) {
    bd=PERIODIC;
}

// bl=(int)(*bleft);
// br=(int)(*bright);
// bu=(int)(*bupper);
// bd=(int)(*blower);

invtwodx=0.5/(*dx);
invtwody=0.5/(*dy);

// interior points
for (j=1; j<(J-1); j++) { /* Loop through columns, excluding x boundaries */

  jdown = j-1; 
  jup = j+1;
  Kj=K*j;
  Kjup=K*jup;
  Kjdown=K*jdown;

  for (i=1; i<(K-1); i++) { /* Loop through rows, excluding y boundaries */              

      gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
      gy=(M[Kj+i+1]-M[Kj+i-1])*invtwody;

      S[Kj+i]=sqrt(gx*gx+gy*gy);

  }
}


// left boundary, except corners
j=0;
Kj=K*j;

switch (bl) {
    case FIXED:
        // slope remains zero
        break;
    case MIRROR: // mirror
        for (i=1; i<(K-1); i++) {
            gx=0;        
            gy=(M[Kj+i+1]-M[Kj+i-1])*invtwody;
            S[Kj+i]=fabs(gy);
        }
        break;
    case PERIODIC: // periodic
        jdown = J-1; 
        jup = j+1;
        Kjup=K*jup;
        Kjdown=K*jdown;
        for (i=1; i<(K-1); i++) {
            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
            gy=(M[Kj+i+1]-M[Kj+i-1])*invtwody;
            S[Kj+i]=sqrt(gx*gx+gy*gy);
        }
        break;
}


// right boundary, except corners
j=J-1;
Kj=K*j;

switch (br) {
    case FIXED:
        // slope remains zero
        break;
    case MIRROR: // mirror
        for (i=1; i<(K-1); i++) {
            gx=0;        
            gy=(M[Kj+i+1]-M[Kj+i-1])*invtwody;
            S[Kj+i]=fabs(gy);
        }
        break;
    case PERIODIC: // periodic
        jdown = j-1; 
        jup = 0;
        Kjup=K*jup;
        Kjdown=K*jdown;
        for (i=1; i<(K-1); i++) {
            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
            gy=(M[Kj+i+1]-M[Kj+i-1])*invtwody;
            S[Kj+i]=sqrt(gx*gx+gy*gy);
        }
        break;
}


// upper boundary, except corners
i=0;

switch (bu) {
    case FIXED:
        // slope remains zero
        break;
    case MIRROR: // mirror
        for (j=1; j<(J-1); j++) {
            jdown = j-1; 
            jup = j+1;
            Kj=K*j;
            Kjup=K*jup;
            Kjdown=K*jdown;

            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;        
            gy=0;
            S[Kj+i]=fabs(gx);
        }
        break;
    case PERIODIC: // periodic
        iup=i+1; idown=K-1;
        for (j=1; j<(J-1); j++) {
            jdown = j-1; jup = j+1;
            Kj=K*j;
            Kjup=K*jup;
            Kjdown=K*jdown;

            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
            gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
            S[Kj+i]=sqrt(gx*gx+gy*gy);
        }
        break;
}


// lower boundary, except corners
i=K-1;

switch (bd) {
    case FIXED:
        // slope remains zero
        break;
    case MIRROR: // mirror
        for (j=1; j<(J-1); j++) {
            jdown = j-1; jup = j+1;
            Kj=K*j;
            Kjup=K*jup;
            Kjdown=K*jdown;

            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;        
            gy=0;
            S[Kj+i]=fabs(gx);
        }
        break;
    case PERIODIC: // periodic
        iup=0; idown=i-1;
        for (j=1; j<(J-1); j++) {
            jdown = j-1; jup = j+1;
            Kj=K*j;
            Kjup=K*jup;
            Kjdown=K*jdown;

            gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
            gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
            S[Kj+i]=sqrt(gx*gx+gy*gy);
        }
        break;
}


// UL corner
i=0;
j=0;

jdown = J-1; jup = j+1;
Kj=K*j;
Kjup=K*jup;
Kjdown=K*jdown;

switch (bl) {
    case FIXED: // fixed x
        // slope remains zero
        break;
    case MIRROR: // mirror x
        gx=0;
        switch (bu) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                // slope remains zero
                break;
            case PERIODIC:
                iup=i+1; idown=K-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=fabs(gy);
                break;
        }
        break;
    case PERIODIC: // periodic x
        gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
        switch (bu) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                gy=0;
                S[Kj+i]=fabs(gx);
                break;
            case PERIODIC:
                iup=i+1; idown=K-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=sqrt(gx*gx+gy*gy);
                break;
        }
        break;
}

// UR corner
i=0;
j=J-1;

jdown = j-1; jup = 0;
Kj=K*j;
Kjup=K*jup;
Kjdown=K*jdown;

switch (br) {
    case FIXED: // fixed x
        // slope remains zero
        break;
    case MIRROR: // mirror x
        gx=0;
        switch (bu) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                // slope remains zero
                break;
            case PERIODIC:
                iup=i+1; idown=K-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=fabs(gy);
                break;
        }
        break;
    case PERIODIC: // periodic x
        gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
        switch (bu) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                gy=0;
                S[Kj+i]=fabs(gx);
                break;
            case PERIODIC:
                iup=i+1; idown=K-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=sqrt(gx*gx+gy*gy);
                break;
        }
        break;
}

// LL corner
i=K-1;
j=0;

jdown = J-1; jup = j+1;
Kj=K*j;
Kjup=K*jup;
Kjdown=K*jdown;

switch (bl) {
    case FIXED: // fixed x
        // slope remains zero
        break;
    case MIRROR: // mirror x
        gx=0;
        switch (bd) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                // slope remains zero
                break;
            case PERIODIC:
                iup=0; idown=i-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=fabs(gy);
                break;
        }
        break;
    case PERIODIC: // periodic x
        gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
        switch (bd) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                gy=0;
                S[Kj+i]=fabs(gx);
                break;
            case PERIODIC:
                iup=0; idown=i-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=sqrt(gx*gx+gy*gy);
                break;
        }
        break;
}


// LR corner
i=K-1;
j=J-1;

jdown = j-1; jup = 0;
Kj=K*j;
Kjup=K*jup;
Kjdown=K*jdown;

switch (br) {
    case FIXED: // fixed x
        // slope remains zero
        break;
    case MIRROR: // mirror x
        gx=0;
        switch (bd) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                // slope remains zero
                break;
            case PERIODIC:
                iup=0; idown=i-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=fabs(gy);
                break;
        }
        break;
    case PERIODIC: // periodic x
        gx=(M[Kjup+i]-M[Kjdown+i])*invtwodx;
        switch (bd) {
            case FIXED:
                // slope remains zero
                break;
            case MIRROR:
                gy=0;
                S[Kj+i]=fabs(gx);
                break;
            case PERIODIC:
                iup=0; idown=i-1;
                gy=(M[Kj+iup]-M[Kj+idown])*invtwody;
                S[Kj+i]=sqrt(gx*gx+gy*gy);
                break;
        }
        break;
}


} /* end GetSlopes subroutine */ 



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *M, *S, *dx, *dy, *b;
    int K, J;
    
    /* Skip argument checking for efficiency */
    
    
    /* Get pointers to inputs */
    M = (double *)mxGetPr(prhs[0]); /* elevations */
    dx = (double *)mxGetPr(prhs[1]); /* dx */
    dy = (double *)mxGetPr(prhs[2]); /* dy */
    b = (double *)mxGetPr(prhs[3]); /* [bleft bright bupper blower] */
    
    /* Get dimensions of input matrix of elevations */
    K = mxGetM(prhs[0]); J = mxGetN(prhs[0]);
    
    /* Create array for the return argument and get a pointer to it */
    S = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); /* matrix of slopes */
    
    GetSlopes(K, J, M, S, dx, dy, b);
    
} /* End mexFunction */
