/*
 * =============================================================
 * mexHillLinear.c - C MEX file for one implicit step of linear
 *                   hillslope diffusion (p.Sc = Inf; HillslopeLinear.m)
 *
 * usage: [Z, Deliv, info] = mexHillLinear(Z, Zguess, dodiff, Kn, par, bvec)
 *
 *    Advances dz/dt = RHS(z) over dt for the cells that diffuse, where RHS is the
 *    conservative 5-point operator of AndrewsHanks / HillslopeFaces with the
 *    diffusivity D = Kn (no slope dependence):
 *
 *      RHS_c = sum over the four faces f of c_f*(z_nbr - z_c),
 *      c_f = (D_c + D_nbr)/(2*h_f^2)
 *
 *    with D = 0 at cells that don't diffuse. On a face between a diffusing
 *    cell and a cell that doesn't diffuse (a sink), the diffusing cell's D is
 *    used instead of the average (c_f = (D_c + D_nbr)/h_f^2), so the flux into
 *    sinks is not halved. Neighbours across mirror edges are the reflected
 *    cells, across periodic edges the wrapped cells, and beyond fixed edges
 *    ghosts of zero elevation and zero D (as in PadGrid). Sinks keep their
 *    elevation.
 *
 *    Backward Euler (method 0), z - dt*RHS(z) = z^n: the unknowns are all
 *    cells; rows of sinks are identity rows,
 *    and their coupling to diffusing cells is moved to the right side. Each
 *    row is multiplied by the fraction of its cell inside the domain (omega:
 *    1/2 on mirror edges, 1/4 at mirror corners), which makes the matrix
 *    symmetric (a mirror-edge cell couples to its interior neighbour across
 *    two faces). The matrix is then a symmetric M-matrix, solved by conjugate
 *    gradients with an incomplete Cholesky (IC(0)) preconditioner on its own
 *    sparsity pattern (periodic links included), starting from Zguess. The
 *    iteration stops when the largest residual per unit cell area is below
 *    tol times the largest change of z, or below 1e-13 of the largest |z^n|.
 *    The new elevation of each diffusing cell is then z^n + dt*RHS(z), and
 *    each sink receives dt*RHS(z), so the volume budget closes to round-off
 *    whatever the tolerance.
 *
 *    Peaceman-Rachford ADI (method 1; Crank-Nicolson in 2D), with RHS split
 *    into its x and y faces, Lx + Ly:
 *      (I - dt/2 Lx) u* = (I + dt/2 Ly) z^n,   (I - dt/2 Ly) u = (I + dt/2 Lx) u*
 *    Each half step is a set of tridiagonal systems along grid lines (cyclic
 *    on periodic lines, by Sherman-Morrison), solved directly (Thomas). O(N),
 *    second order in time, but not L-stable: at dt*D/dx^2 >> 1 it can overshoot
 *    next to sinks. The new elevation of each
 *    diffusing cell is z^n + dt*Lx u* + dt/2*Ly (z^n + u) from the face fluxes,
 *    and sinks receive the same, so the budget closes to round-off.
 *
 *    Input arguments (K x J matrices unless noted):
 *      Z:       elevation at the start of the step
 *      Zguess:  starting guess for the solution ([] to start from Z)
 *      dodiff:  1 where cells diffuse, 0 at sinks (fixed points), 2 at channel
 *               cells that don't diffuse: sinks whose faces to lower diffusing
 *               neighbours are closed (backward Euler only; the solve is repeated
 *               until those faces agree with the solution; hilldiffusion.h) (double)
 *      Kn:      diffusivity, scalar or K x J
 *      par:     [dx dy dt tol maxIter method]; method 0 (default): backward
 *               Euler; 1: Peaceman-Rachford ADI (tol, maxIter and Zguess unused)
 *      bvec:    boundary codes [left right upper lower] (0 fixed, 1 mirror, 2 periodic)
 *
 *    Return arguments:
 *      Z:       elevation at the end of the step
 *      Deliv:   thickness of material delivered to sinks during the step (0 elsewhere)
 *      info:    [iterations converged maxResidual resolves] ([0 1 0 0] with ADI;
 *               resolves: extra solves for the faces of dodiff = 2 cells)
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2026 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include "hilldiffusion.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    const double *Z, *Zg = NULL, *dodiff, *Kn, *par, *bvec;
    double *Zn, *Deliv, *info, *cf, *om, dx, dy, dt, tol, maxF = 0;
    int K, J, N, i, c, sKn, maxIter, method = 0, converged = 1, iter = 0, resolves = 0;
    int *nbr;

    if (nrhs != 6)
        mexErrMsgIdAndTxt("mexHillLinear:nrhs", "Usage: [Z, Deliv, info] = mexHillLinear(Z, Zguess, dodiff, Kn, par, bvec)");
    for (i = 0; i < nrhs; i++)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
            mexErrMsgIdAndTxt("mexHillLinear:type", "Input %d must be a real double array.", i + 1);
    if (mxGetNumberOfElements(prhs[4]) < 5)
        mexErrMsgIdAndTxt("mexHillLinear:par", "par must be [dx dy dt tol maxIter].");

    Z      = mxGetPr(prhs[0]);
    dodiff = mxGetPr(prhs[2]);
    Kn     = mxGetPr(prhs[3]);
    par    = mxGetPr(prhs[4]);
    bvec   = mxGetPr(prhs[5]);
    K = (int)mxGetM(prhs[0]);
    J = (int)mxGetN(prhs[0]);
    N = K*J;
    if (mxGetNumberOfElements(prhs[1]) == (size_t)N) Zg = mxGetPr(prhs[1]);
    if ((int)mxGetNumberOfElements(prhs[2]) != N)
        mexErrMsgIdAndTxt("mexHillLinear:size", "dodiff must be K x J.");
    sKn = mxGetNumberOfElements(prhs[3]) > 1; /* stride: 0 for a scalar, 1 for a matrix */

    dx = par[0]; dy = par[1]; dt = par[2]; tol = par[3]; maxIter = (int)par[4];
    if (mxGetNumberOfElements(prhs[4]) > 5) method = (int)par[5];

    Zn    = mxGetPr(plhs[0] = mxCreateDoubleMatrix(K, J, mxREAL));
    Deliv = mxGetPr(plhs[1] = mxCreateDoubleMatrix(K, J, mxREAL));
    info  = mxGetPr(plhs[2] = mxCreateDoubleMatrix(1, 4, mxREAL));

    /* face neighbours, face conductances, and the fraction of each cell inside the domain */
    nbr = (int *)mxMalloc(4*N*sizeof(int));
    cf  = (double *)mxMalloc(4*N*sizeof(double));
    om  = (double *)mxMalloc(N*sizeof(double));
    HDFaces(K, J, dodiff, Kn, sKn, dx, dy, bvec, nbr, cf, om);

    if (method == 1) {
        HDSolveADI(K, J, Z, dodiff, nbr, cf, bvec, dt, Zn, Deliv);
    } else {
        HDWork w;
        int has2 = 0;
        HDAlloc(&w, N);
        iter = HDSolveBE(&w, N, Z, Zg, NULL, dodiff, nbr, cf, om, dt, tol, maxIter, Zn, Deliv, &converged, &maxF);
        for (c = 0; c < N; c++) if (dodiff[c] == 2) { has2 = 1; break; }
        if (has2) {
            /* channel cells (dodiff = 2; hilldiffusion.h): the closed faces were
               decided on Zg (or z^n); solve again with the faces of the solution,
               starting from it, until they agree */
            double *zref = (double *)mxMalloc(N*sizeof(double));
            const double *zr = Zg ? Zg : Z;
            int k;
            for (k = 0; k < 20 && HDClosedChanged(N, nbr, dodiff, Z, zr, Zn) > 0; k++) {
                memcpy(zref, Zn, N*sizeof(double));
                zr = zref;
                iter += HDSolveBE(&w, N, Z, zref, zref, dodiff, nbr, cf, om, dt, tol, maxIter, Zn, Deliv, &converged, &maxF);
                resolves++;
            }
            if (k == 20 && HDClosedChanged(N, nbr, dodiff, Z, zr, Zn) > 0)
                mexWarnMsgIdAndTxt("Tadpole:hillLinear", "Faces of channel cells did not settle in 20 solves.");
            mxFree(zref);
        }
        HDFree(&w);
    }
    info[0] = iter;
    info[1] = converged;
    info[2] = maxF;
    info[3] = resolves; /* extra solves for the faces of channel cells */

    mxFree(nbr); mxFree(cf); mxFree(om);
}
