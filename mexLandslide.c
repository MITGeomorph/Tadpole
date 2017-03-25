 /*
 * =============================================================
 * mexLandslide.c - C MEX file to reduce slopes to angle of repose
 *
 * usage: [M dz] = mexLandslide(M,order,dx,dy,Sc,bvec,fixed)
 *
 *    where M is a K x J matrix of elevations and dx, dy are the cell
 *    dimensions. 
 *    bvec is a vector of codes for the [left, right, upper, lower]
 *    boundary conditions: 0=fixed, 1=mirror, 2=periodic
 *
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


// Landslide(o[i]-1,K,J,M,dzmax,b,f,dz);
void Landslide(double ind, const int K, const int J, double M[], double dzmax[], double bdy[], double f[], double dz[])
{

int i, j, k, p, q, n, bl, br, bu, bd;
int di[] = {0,-1,-1,-1,0,1,1,1};
int dj[] = {1,1,0,-1,-1,-1,0,1};
int skip[] = {0,0,0,0,0,0,0,0};
double e0, en;

// convert linear index to row,col
n = (int)ind;

i = n % K;
j = (n-i)/K;



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
        f[n] = 1;
        break;
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
        skip[5]=1;
        skip[6]=1;
        skip[7]=1;
        f[n] = 1;
        break;
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
        skip[3]=1;
        skip[4]=1;
        skip[5]=1;
        f[n] = 1;
        break;
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
        skip[7]=1;
        skip[0]=1;
        skip[1]=1;
        f[n] = 1;
        break;
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


if (!f[n]) { // leave fixed points alone. Note that we made sure this includes fixed boundaries in the switch statements above
    
    e0 = M[n]; // elevation at i,j
    
//     minima[n] = 1; // minimum until proven otherwise
    
    for (k=0; k<8; k++) {  /* loop through neighbors */

        if (!skip[k]) {
            
            p = i + di[k];
            q = j + dj[k];
            
            en = M[K*q+p]; // elevation of neighbor k
                        
//             if (e0 > en) { // if it's a downslope neighbor
//                 minima[n] = 0; // there is at least one downslope neighbor, so it's not a minimum
                if ((e0 - en) > dzmax[k]) { // if the slope is steeper than Sc
                    e0 = en + dzmax[k]; // set it to Sc
                }
//             }
        }
        
    } /* end of neighbors loop */
    
    // record the elevation change and the new elevation
    dz[n] = M[n] - e0;
    M[n] = e0;
}

} /* end Landslide subfunction */




//  usage: [M dz] = mexLandslide(M,order,dx,dy,Sc,bvec,fixed)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *Min, *Mout, *dz, *o, *dx, *dy, *Sc, *b, *f, dzmax[8];
//   double *minima;
  int i, k, K, J;

  // Get pointers to inputs
  Min = (double *)mxGetPr(prhs[0]); // elevations
  o = (double *)mxGetPr(prhs[1]); // array giving the one-based ordering of the elements in M, in order of increasing elevation
  dx = (double *)mxGetPr(prhs[2]); 
  dy = (double *)mxGetPr(prhs[3]); 
  Sc = (double *)mxGetPr(prhs[4]); // angle of repose (m/m)
  b = (double *)mxGetPr(prhs[5]); /* [bleft bright bupper blower] */
  f = (double *)mxGetPr(prhs[6]); /* binary matrix indicating fixed points */
  
  // Get dimensions of input matrix of elevations
  K = mxGetM(prhs[0]);
  J = mxGetN(prhs[0]);

  // pointers to outputs
  Mout = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); // the post-landslide elevations
//   minima = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // locations of local minima
  dz = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); // points that changed elevation

  // copy the input elevations. We do this to avoid modifying the input in memory (efficient, but risky)
  for (i=0; i<(K*J); i++) {
      Mout[i] = Min[i];
  }
  
  // calculate maximum slope in each direction
  for (k=0; k<8; k++) {
      switch (k) {
          case 0:
          case 4:
              dzmax[k] = (*Sc)*(*dx);
              break;
          case 2:
          case 6:
              dzmax[k] = (*Sc)*(*dy);
              break;
          case 1:
          case 3:
          case 5:
          case 7:
              dzmax[k] = (*Sc)*sqrt((*dx)*(*dx) + (*dy)*(*dy));
              break;
      }
  }
  
  // adjust elevation of each element, in order of increasing elevation
  for (i=0; i<(K*J); i++) { // all points
//       Landslide(o[i]-1,K,J,Mout,dzmax,b,f,minima);
      Landslide(o[i]-1,K,J,Mout,dzmax,b,f,dz);
  }

} // End mexFunction
