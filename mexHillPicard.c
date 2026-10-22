/*
 * =============================================================
 * mexHillPicard.c - C MEX file for one fully implicit (backward Euler) step of
 *                   nonlinear hillslope diffusion, solved by Picard iteration
 *                   with Anderson acceleration
 *
 * usage: [Z, Deliv, info, hist] = mexHillPicard(Z, Zguess, dodiff, Knl, par, bvec)
 *
 *    Solves z - dt*RHS(z) = z^n, where RHS is the conservative 5-point operator
 *    of AndrewsHanks / HillslopeFaces with the cell diffusivity
 *
 *      D = Knl*(1 + s + s^2 + ... + s^N),   s = (S/Sc)^2,
 *
 *    S from centred differences (padded as in PadGrid), capped at cap*max(Knl)
 *    (the series approximation of Ganti et al., 2012, as in the explicit
 *    solver; with N = 5 and cap = 6, the value at Sc, it is the explicit solver's law). Picard
 *    iteration: z^(k+1) = solution of z - dt*L(D(z^k)) z = z^n, each a linear
 *    backward Euler solve (hilldiffusion.h: IC(0)-preconditioned conjugate
 *    gradients, warm-started from the current iterate). Plain Picard can cycle
 *    where slopes are near or above Sc (D jumps between large and small
 *    values at the same cells), so the iterates are combined by Anderson
 *    acceleration over the last m iterates: with f_k = G(z_k) - z_k,
 *    z_(k+1) = G(z_k) - dG*gamma, gamma minimizing |f_k - dF*gamma| (dF, dG
 *    the differences of successive f and G), solved by the normal equations.
 *    The iteration stops when the largest change of z between iterations is
 *    below tol times the largest change in the step. The result is the last
 *    linear solve, whose new elevations and delivery to sinks come from the
 *    face fluxes, so volume is conserved to round-off whatever the tolerance.
 *
 *    Input arguments (K x J matrices unless noted):
 *      Z:       elevation at the start of the step
 *      Zguess:  starting guess ([] to start from Z)
 *      dodiff:  1 where cells diffuse, 0 at sinks (double)
 *      Knl:     hillslope diffusivity, scalar or K x J
 *      par:     [dx dy dt Sc N cap tol maxIter m cgTol]
 *      bvec:    boundary codes [left right upper lower] (0 fixed, 1 mirror, 2 periodic)
 *
 *    Return arguments:
 *      Z:       elevation at the end of the step
 *      Deliv:   thickness of material delivered to sinks (0 elsewhere)
 *      info:    [iterations converged cgIterations lastChange]
 *      hist:    (optional) largest change in each iteration (1 x iterations)
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2026 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include "hilldiffusion.h"


/* cell diffusivity: Knl*(1 + s + ... + s^N), capped at Dcap */
static void CellD(int K, int J, const double *z, const double *Knl, int sK, double dx, double dy,
                  double Sc, int N, double Dcap, const double *bvec, double *D) {
    int i, j, c, k;
    for (j = 0; j < J; j++)
        for (i = 0; i < K; i++) {
            int nr = FaceNeighbor(i, j, 0, K, J, bvec), nl = FaceNeighbor(i, j, 1, K, J, bvec);
            int nu = FaceNeighbor(i, j, 2, K, J, bvec), nd = FaceNeighbor(i, j, 3, K, J, bvec);
            double zx, zy, s, term, series;
            c = i + K*j;
            zx = ((nr >= 0 ? z[nr] : 0) - (nl >= 0 ? z[nl] : 0))/(2*dx);
            zy = ((nu >= 0 ? z[nu] : 0) - (nd >= 0 ? z[nd] : 0))/(2*dy);
            s = (zx*zx + zy*zy)/(Sc*Sc);
            term = 1; series = 1;
            for (k = 1; k <= N; k++) { term *= s; series += term; }
            D[c] = Knl[sK*c]*series;
            if (D[c] > Dcap) D[c] = Dcap;
        }
}


/* Solve the m x m system M*g = v (symmetric positive semi-definite) by
   Cholesky, with a small diagonal shift for nearly dependent columns. */
static void SolveSmall(int m, double *M, double *v, double *g) {
    double L[100], shift = 0, y[10];
    int i, j, k;
    for (i = 0; i < m; i++) if (M[i*m+i] > shift) shift = M[i*m+i];
    shift *= 1e-12;
    for (i = 0; i < m; i++) {
        for (j = 0; j <= i; j++) {
            double s = M[i*m+j] + (i == j ? shift : 0);
            for (k = 0; k < j; k++) s -= L[i*m+k]*L[j*m+k];
            if (i == j) L[i*m+i] = sqrt(s > 0 ? s : 1e-300);
            else L[i*m+j] = s/L[j*m+j];
        }
    }
    for (i = 0; i < m; i++) {
        double s = v[i];
        for (k = 0; k < i; k++) s -= L[i*m+k]*y[k];
        y[i] = s/L[i*m+i];
    }
    for (i = m-1; i >= 0; i--) {
        double s = y[i];
        for (k = i+1; k < m; k++) s -= L[k*m+i]*g[k];
        g[i] = s/L[i*m+i];
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    const double *Z, *Zg = NULL, *dodiff, *Knl, *par, *bvec;
    double *Zn, *Deliv, *info, *hist, *cf, *om, *D, *z, *dF, *dG, *fprev, *gprev;
    double dx, dy, dt, Sc, cap, tol, cgTol, Dcap, Kmax = 0, zscale = 1, change = 0;
    int K, J, N, i, c, sK, Nser, maxIter, m, it, converged = 0, cgTotal = 0, nd, ncol = 0, head = 0;
    int *nbr, *dlist;
    HDWork w;

    if (nrhs != 6)
        mexErrMsgIdAndTxt("mexHillPicard:nrhs", "Usage: [Z, Deliv, info, hist] = mexHillPicard(Z, Zguess, dodiff, Knl, par, bvec)");
    for (i = 0; i < nrhs; i++)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
            mexErrMsgIdAndTxt("mexHillPicard:type", "Input %d must be a real double array.", i + 1);
    if (mxGetNumberOfElements(prhs[4]) < 10)
        mexErrMsgIdAndTxt("mexHillPicard:par", "par must be [dx dy dt Sc N cap tol maxIter m cgTol].");

    Z      = mxGetPr(prhs[0]);
    dodiff = mxGetPr(prhs[2]);
    Knl    = mxGetPr(prhs[3]);
    par    = mxGetPr(prhs[4]);
    bvec   = mxGetPr(prhs[5]);
    K = (int)mxGetM(prhs[0]);
    J = (int)mxGetN(prhs[0]);
    N = K*J;
    if (mxGetNumberOfElements(prhs[1]) == (size_t)N) Zg = mxGetPr(prhs[1]);
    if ((int)mxGetNumberOfElements(prhs[2]) != N)
        mexErrMsgIdAndTxt("mexHillPicard:size", "dodiff must be K x J.");
    sK = mxGetNumberOfElements(prhs[3]) > 1;

    dx = par[0]; dy = par[1]; dt = par[2]; Sc = par[3]; Nser = (int)par[4]; cap = par[5];
    tol = par[6]; maxIter = (int)par[7]; m = (int)par[8]; cgTol = par[9];
    if (m > 10) m = 10;
    for (c = 0; c < N; c++) {
        if (Knl[sK*c] > Kmax) Kmax = Knl[sK*c];
        if (fabs(Z[c]) > zscale) zscale = fabs(Z[c]);
    }
    Dcap = cap*Kmax;

    Zn    = mxGetPr(plhs[0] = mxCreateDoubleMatrix(K, J, mxREAL));
    Deliv = mxGetPr(plhs[1] = mxCreateDoubleMatrix(K, J, mxREAL));
    info  = mxGetPr(plhs[2] = mxCreateDoubleMatrix(1, 4, mxREAL));
    hist  = (double *)mxCalloc(maxIter + 1, sizeof(double));

    nbr = (int *)mxMalloc(4*N*sizeof(int));
    cf  = (double *)mxMalloc(4*N*sizeof(double));
    om  = (double *)mxMalloc(N*sizeof(double));
    D   = (double *)mxMalloc(N*sizeof(double));
    z   = (double *)mxMalloc(N*sizeof(double));
    HDAlloc(&w, N);

    /* diffusing cells: the unknowns of the Anderson combination */
    dlist = (int *)mxMalloc(N*sizeof(int));
    nd = 0;
    for (c = 0; c < N; c++) if (dodiff[c] == 1) dlist[nd++] = c; /* channel cells (dodiff = 2) keep their elevation */
    dF = (double *)mxMalloc(((size_t)nd*(m > 0 ? m : 1) + 1)*sizeof(double));
    dG = (double *)mxMalloc(((size_t)nd*(m > 0 ? m : 1) + 1)*sizeof(double));
    fprev = (double *)mxMalloc((nd + 1)*sizeof(double));
    gprev = (double *)mxMalloc((nd + 1)*sizeof(double));

    for (c = 0; c < N; c++) z[c] = (Zg && dodiff[c] == 1) ? Zg[c] : Z[c];

    for (it = 1; it <= maxIter; it++) {
        int cgc, k, l;
        double maxF, dzmax = 0;

        /* G(z): linear backward Euler solve with D(z) */
        CellD(K, J, z, Knl, sK, dx, dy, Sc, Nser, Dcap, bvec, D);
        HDFaces(K, J, dodiff, D, 1, dx, dy, bvec, nbr, cf, om);
        cgTotal += HDSolveBE(&w, N, Z, z, z, dodiff, nbr, cf, om, dt, cgTol, 1000, Zn, Deliv, &cgc, &maxF);

        change = 0;
        for (c = 0; c < N; c++) {
            double d = fabs(Zn[c] - z[c]), dz = fabs(Zn[c] - Z[c]);
            if (d > change) change = d;
            if (dz > dzmax) dzmax = dz;
        }
        hist[it-1] = change;
        if (change <= tol*dzmax || change <= 1e-12*zscale) { converged = 1; break; }

        if (m == 0) {
            memcpy(z, Zn, N*sizeof(double));
            continue;
        }

        /* Anderson: f = G(z) - z, g = G(z); differences of successive f and g
           kept in a circular buffer of m columns */
        if (it > 1) {
            double *colF = dF + (size_t)head*nd, *colG = dG + (size_t)head*nd;
            for (k = 0; k < nd; k++) {
                c = dlist[k];
                colF[k] = (Zn[c] - z[c]) - fprev[k];
                colG[k] = Zn[c] - gprev[k];
            }
            head = (head + 1) % m;
            if (ncol < m) ncol++;
        }
        for (k = 0; k < nd; k++) {
            c = dlist[k];
            fprev[k] = Zn[c] - z[c];
            gprev[k] = Zn[c];
        }
        if (ncol == 0) {
            memcpy(z, Zn, N*sizeof(double));
        } else {
            double M[100], v[10], gam[10];
            for (k = 0; k < ncol; k++) {
                const double *Fk = dF + (size_t)k*nd;
                v[k] = 0;
                for (i = 0; i < nd; i++) v[k] += Fk[i]*fprev[i];
                for (l = 0; l <= k; l++) {
                    const double *Fl = dF + (size_t)l*nd;
                    double s = 0;
                    for (i = 0; i < nd; i++) s += Fk[i]*Fl[i];
                    M[k*ncol+l] = s; M[l*ncol+k] = s;
                }
            }
            SolveSmall(ncol, M, v, gam);
            memcpy(z, Zn, N*sizeof(double)); /* sinks: Z */
            for (i = 0; i < nd; i++) {
                double s = gprev[i];
                for (k = 0; k < ncol; k++) s -= dG[(size_t)k*nd + i]*gam[k];
                z[dlist[i]] = s;
            }
        }
    }
    if (it > maxIter) it = maxIter;
    if (!converged)
        mexWarnMsgIdAndTxt("Tadpole:hillPicard",
            "Picard iteration did not converge in %d iterations (max change %g m).", maxIter, change);

    info[0] = it;
    info[1] = converged;
    info[2] = cgTotal;
    info[3] = change;
    if (nlhs > 3) {
        double *h = mxGetPr(plhs[3] = mxCreateDoubleMatrix(1, it, mxREAL));
        memcpy(h, hist, it*sizeof(double));
    }

    HDFree(&w);
    mxFree(nbr); mxFree(cf); mxFree(om); mxFree(D); mxFree(z); mxFree(dlist);
    mxFree(dF); mxFree(dG); mxFree(fprev); mxFree(gprev); mxFree(hist);
}
