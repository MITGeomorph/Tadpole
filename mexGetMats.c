// mexGetMats.c
// [Ar Ac Av Br Bc Bv amax bmax] = mexGetMats(D,dt,dx,dy,bdy,C)

#include "mex.h"
#include "matrix.h"

#define FIXED    0
#define MIRROR   1
#define PERIODIC 2


// nnz = BuildMats(K,J,Ar,Ac,Av,Br,Bc,Bv,D,dt,dx,dy,bdy,C,amax,bmax)
void BuildMats(const int K, const int J, double Ar[], double Ac[], double Av[], double Br[], double Bc[], double Bv[], double *D, double *dt, double *dx, double *dy, double bdy[], double C[], double *amax, double *bmax)
{

int i, j, a, b, c, N = K*J;
int bl, br, bu, bd;
double Xn = (*D)*(*dt)/(2*(*dx)*(*dx)); 
double Xc = -2*Xn;
double Yn = (*D)*(*dt)/(2*(*dy)*(*dy)); 
double Yc = -2*Yn;

        
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

    

a = 0; // keeps track of which element we're at in the A vectors
b = 0; // keeps track of which element we're at in the B vectors


// to do: 
// x change loops
// x change 1-based to zero-based
// x array indexing with [c] rather than (c)
// x switch statements, including break;
// x a++

// Interior points
for (i=2; i<K; i++) { // loop through interior elements
   for (j=2; j<J; j++) {
       c = K*(j-1)+i; // index of current element
       if (C[c-1]) {

          // diagonal elements
          Ar[a]=c;
          Ac[a]=c;
          Av[a]=Xc;
          a++;

          Br[b]=c;
          Bc[b]=c;
          Bv[b]=Yc;
          b++;

          // off-diagonal elements

          // W neighbor
          Ar[a]=c;
          Ac[a]=c-K;
          Av[a]=Xn;
          a++;

          // E neighbor
          Ar[a]=c;
          Ac[a]=c+K;
          Av[a]=Xn;
          a++;

          // N neighbor
          Br[b]=c;
          Bc[b]=c-1;
          Bv[b]=Yn;
          b++;

          // S neighbor
          Br[b]=c;
          Bc[b]=c+1;
          Bv[b]=Yn;
          b++;
       }
   }
}

// x derivatives on boundaries:

// left boundary
j=1;
switch (bl) {

    case FIXED:
        // leave zeros
        break;
    case MIRROR:
        for (i=1; i<(K+1); i++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal element
                Ar[a]=c;
                Ac[a]=c;
                Av[a]=Xc;
                a++;

                // left neighbor is off the grid

                // right neighbor
                Ar[a]=c;
                Ac[a]=c+K;
                Av[a]=2*Xn; // the mirror condition
                a++;

            }
        }
        break;
    case PERIODIC:
        for (i=1; i<(K+1); i++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal element
                Ar[a]=c;
                Ac[a]=c;
                Av[a]=Xc;
                a++;

                // left neighbor
                Ar[a]=c;
                Ac[a]=c-K+N; // the periodic condition
                Av[a]=Xn;
                a++;

                // right neighbor
                Ar[a]=c;
                Ac[a]=c+K;
                Av[a]=Xn;
                a++;
            }
        }
        break;
}

// right boundary
j=J;
switch (br) {

    case FIXED:
        // leave zeros
        break;
    case MIRROR:
        for (i=1; i<(K+1); i++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal element
                Ar[a]=c;
                Ac[a]=c;
                Av[a]=Xc;
                a++;

                // left neighbor
                Ar[a]=c;
                Ac[a]=c-K; // the periodic condition
                Av[a]=2*Xn;
                a++;

                // right neighbor is off the grid
            }
        }
        break;
    case PERIODIC:
        for (i=1; i<(K+1); i++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal element
                Ar[a]=c;
                Ac[a]=c;
                Av[a]=Xc;
                a++;

                // left neighbor
                Ar[a]=c;
                Ac[a]=c-K; 
                Av[a]=Xn;
                a++;

                // right neighbor
                Ar[a]=c;
                Ac[a]=c+K-N; // the periodic condition
                Av[a]=Xn;
                a++;
            }
        }
        break;
}


// Upper boundary excluding corners
i=1;
if (bu == FIXED) {
    // Leave zeros 
} else {
   for (j=2; j<J; j++) {
        c=K*(j-1)+i; // index of current element
        if (C[c-1]) {

          // diagonal elements
          Ar[a]=c;
          Ac[a]=c;
          Av[a]=Xc;
          a++;

          // off-diagonal elements

          // left neighbor
          Ar[a]=c;
          Ac[a]=c-K;
          Av[a]=Xn;
          a++;

          // right neighbor
          Ar[a]=c;
          Ac[a]=c+K;
          Av[a]=Xn;
          a++;
        }       
   }
}

// Lower boundary excluding corners
i=K;
if (bd == FIXED) {
    // leave zero 
} else {
   for (j=2; j<J; j++) {
        c=K*(j-1)+i; // index of current element
        if (C[c-1]) {

          // diagonal elements
          Ar[a]=c;
          Ac[a]=c;
          Av[a]=Xc;
          a++;

          // off-diagonal elements

          // left neighbor
          Ar[a]=c;
          Ac[a]=c-K;
          Av[a]=Xn;
          a++;

          // right neighbor
          Ar[a]=c;
          Ac[a]=c+K;
          Av[a]=Xn;
          a++;
        }       
   }
}




// y derivatives on boundaries:

// upper boundary
i=1;
switch (bu) {

    case FIXED:
        // leave zeros
        break;
    case MIRROR:
        for (j=1; j<(J+1); j++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal elements      
                Br[b]=c;
                Bc[b]=c;
                Bv[b]=Yc;
                b++;

                // off-diagonal elements

                // upper neighbor is off the grid

                // lower neighbor
                Br[b]=c;
                Bc[b]=c+1;
                Bv[b]=2*Yn; // the mirror condition
                b++;
            }
        }
        break;
    case PERIODIC:
        for (j=1; j<(J+1); j++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal elements      
                Br[b]=c;
                Bc[b]=c;
                Bv[b]=Yc;
                b++;

                // off-diagonal elements

                // upper neighbor
                Br[b]=c;
                Bc[b]=c-1+K; // the periodic condition
                Bv[b]=Yn;
                b++;

                // lower neighbor
                Br[b]=c;
                Bc[b]=c+1;
                Bv[b]=Yn;
                b++;
            }
        }
        break;
}

// lower boundary
i=K;
switch (bd) {

    case FIXED:
        // leave zeros
        break;
    case MIRROR:
        for (j=1; j<(J+1); j++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal elements      
                Br[b]=c;
                Bc[b]=c;
                Bv[b]=Yc;
                b++;

                // off-diagonal elements

                // upper neighbor
                Br[b]=c;
                Bc[b]=c-1;
                Bv[b]=2*Yn; // the mirror condition
                b++;

                // lower neighbor is off the grid
            }
        }
        break;
    case PERIODIC:
        for (j=1; j<(J+1); j++) {
            c=K*(j-1)+i; // index of current element
            if (C[c-1]) {

                // diagonal elements      
                Br[b]=c;
                Bc[b]=c;
                Bv[b]=Yc;
                b++;

                // off-diagonal elements

                // upper neighbor
                Br[b]=c;
                Bc[b]=c-1;
                Bv[b]=Yn;
                b++;

                // lower neighbor
                Br[b]=c;
                Bc[b]=c+1-K; // the periodic condition
                Bv[b]=Yn;
                b++;
            }
        }
        break;
}



// left boundary excluding corners
j=1;
if (bl == FIXED) { 
    // Leave zero
} else {
   for (i=2; i<K; i++) { 
        c=K*(j-1)+i; // index of current element
        if (C[c-1]) {

            // diagonal elements      
            Br[b]=c;
            Bc[b]=c;
            Bv[b]=Yc;
            b++;

            // off-diagonal elements

            // upper neighbor
            Br[b]=c;
            Bc[b]=c-1;
            Bv[b]=Yn;
            b++;

            // lower neighbor
            Br[b]=c;
            Bc[b]=c+1;
            Bv[b]=Yn;
            b++;
        }   
   }
    
}

// right boundary excluding corners
j=J;
if (br == FIXED) {
    // Leave zero
} else {
   for (i=2; i<K; i++) { 
        c=K*(j-1)+i; // index of current element
        if (C[c-1]) {

            // diagonal elements      
            Br[b]=c;
            Bc[b]=c;
            Bv[b]=Yc;
            b++;

            // off-diagonal elements

            // upper neighbor
            Br[b]=c;
            Bc[b]=c-1;
            Bv[b]=Yn;
            b++;

            // lower neighbor
            Br[b]=c;
            Bc[b]=c+1;
            Bv[b]=Yn;
            b++;
        }   
   }
}

// the last nonzero elements in 1-based indices
*amax = (double)a; 
*bmax = (double)b;

// eliminate unused elements -- may have to do this back in Matlab
// Aused = (Av~=0);
// Ar = Ar(Aused);
// Ac = Ac(Aused);
// Av = Av(Aused);
// 
// Bused = (Bv~=0);
// Br = Br(Bused);
// Bc = Bc(Bused);
// Bv = Bv(Bused);
        
    
} // end BuildMats


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//   mxArray *fptr, *Fptr;
//   double *z,*bdy,*dt,*dx,*dy,*Nx,*Ny,*a,*b,*U,*f,*F,*Ar,*Ac,*Av,*Af,*B;
//   int K, J, nz;
  double *D, *dt, *dx, *dy, *bdy, *C, *Ar, *Ac, *Av, *Br, *Bc, *Bv, *amax, *bmax;  
  int K, J, nnz;
  
  // Get pointers to inputs
  // syntax is [Ar Ac Av Br Bc Bv] = GetMats(D,dt,dx,dy,bdy,C)
  D = (double *)mxGetPr(prhs[0]); 
  dt = (double *)mxGetPr(prhs[1]);   
  dx = (double *)mxGetPr(prhs[2]); 
  dy = (double *)mxGetPr(prhs[3]); 
  bdy = (double *)mxGetPr(prhs[4]); // 0 fixed, 1 mirror, 2 periodic
  C = (double *)mxGetPr(prhs[5]); // 1 where elevation is allowed to change, 0 not
//   Nx = (double *)mxGetPr(prhs[5]); 
//   Ny = (double *)mxGetPr(prhs[6]); 
//   a = (double *)mxGetPr(prhs[7]); // K/(rho_r/rho_s)
//   b = (double *)mxGetPr(prhs[8]); // 1/Sc^2
//   U = (double *)mxGetPr(prhs[9]);

//   // convert dimensions to integers
//   K = (int)(*Ny+0.5);
//   J = (int)(*Nx+0.5);

  amax = mxCalloc(1, sizeof(int)); // allocate memory for this 1-element integer array
  bmax = mxCalloc(1, sizeof(int)); // allocate memory for this 1-element integer array

  
  K = mxGetM(prhs[5]);
  J = mxGetN(prhs[5]);

  
  // Create arrays for return arguments
  nnz = 3; // Max possible number of nonzero elements per row of the arrays
  Ar = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // 1-BASED (Matlab) row indices 
  Ac = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // 1-BASED (Matlab) column indices 
  Av = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // values  
  Br = (double *)mxGetPr(plhs[3]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // 1-BASED (Matlab) row indices 
  Bc = (double *)mxGetPr(plhs[4]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // 1-BASED (Matlab) column indices
  Bv = (double *)mxGetPr(plhs[5]= mxCreateDoubleMatrix(K*J*nnz, 1, mxREAL)); // values 
  amax = (double *)mxGetPr(plhs[6]= mxCreateDoubleMatrix(1, 1, mxREAL)); // index of last nonzero element in the A vectors 
  bmax = (double *)mxGetPr(plhs[7]= mxCreateDoubleMatrix(1, 1, mxREAL)); // index of last nonzero element in the B vectors 

//   // Create internally used arrays
//   f = (double *)mxGetPr(fptr= mxCreateDoubleMatrix(K*J, 1, mxREAL)); // RHS of time derivative at present time
//   F = (double *)mxGetPr(Fptr= mxCreateDoubleMatrix(K*J, 9, mxREAL)); // partial derivatives of RHS of time derivative at present time
// 
//   GetfF(z,bdy,dx,dy,K,J,a,b,U,f,F); // get values of finite difference approximation and its partial derivatives
  
  BuildMats(K,J,Ar,Ac,Av,Br,Bc,Bv,D,dt,dx,dy,bdy,C,amax,bmax); // calculate values of the nonzero elements in the matrix equation
  
//   *Af = nnz+1; // now *Af is the 1-BASED (Matlab) index of the last nonzero element that was assigned to the A vectors
    
//   // free memory
//   mxDestroyArray(fptr);
//   mxDestroyArray(Fptr);
  
} // End mexFunction
