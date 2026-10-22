/*
 * =============================================================
 * mexSPACE.c - C MEX file for one implicit time step of river incision
 *              and, optionally, sediment entrainment and deposition
 *
 * usage: [Z, R, H, Qs, QsIn, info, lakeFrac] = mexSPACE(Z, R, H, QsPrev, W, Stack,
 *                                  Q, Zf, Lake, erodes, C, Kr, Ks, Vs, par, bvec, QsLat)
 *
 *    Physics (SPACE; Shobe et al., 2017), with Omega = Q^m S^n:
 *      bedrock erosion       Er = Kr Omega exp(-H/H*)        (bed lowering rate)
 *      sediment entrainment  Es = Ks Omega (1 - exp(-H/H*))  (bed lowering rate)
 *      deposition            D  = Vs Qs/Q                    (solid volume per unit bed area and time)
 *      dR/dt = -Er,  dH/dt = D/(1-phi) - Es                  (uplift is applied in Source)
 *    Sediment fluxes Qs are solid volumes. A fraction Ff of eroded bedrock
 *    becomes wash load and leaves the model.
 *
 *    Numerics: flow routing weights are held fixed
 *    during the step. Each cell's new elevation is found implicitly from its
 *    receivers' new elevations, visiting cells receivers first (backward
 *    Euler; Braun & Willett, 2013). The deposition term depends on the
 *    sediment flux from upstream, which is passed downstream in a second
 *    sweep; the two sweeps are repeated until the elevations converge
 *    (Gauss-Seidel iteration, Yuan et al., 2019). Where deposition is strong
 *    (large Vs/P) and the time step long, the inflow and the elevations of a
 *    cell and its donors can feed back on each other so that the iteration
 *    oscillates; the inflow used for deposition is then under-relaxed between
 *    iterations (see par below). The flux leaving each cell
 *    is computed from its actual changes in R and H, so sediment mass is
 *    conserved exactly whatever the convergence tolerance. Sediment entering
 *    a lake raises every lake cell by the same fraction of its water depth.
 *    Without sediment, a single receivers-first sweep is made (implicit
 *    stream power) and only Z is changed.
 *
 *    Input arguments (K x J matrices unless noted):
 *      Z, R, H:  elevation, bedrock elevation and sediment thickness at the
 *                start of the step (R and H may be [] without sediment)
 *      QsPrev:   sediment flux leaving each cell in the previous step; routed
 *                through the current weights, it is the starting guess for the
 *                inflow to each cell (may be [] without sediment)
 *      W:        K x J x 8 flow routing weights (from mexD8, mexDinf, mexDms)
 *      Stack:    (K*J) x 1 1-based linear indices, each cell before its donors
 *      Q:        water discharge, Pav.*A (m^3/yr)
 *      Zf:       filled surface (lake water level)
 *      Lake:     lake label of each flooded cell (its spill cell), 0 elsewhere
 *      erodes:   1 where rivers may erode (channels that are not flooded and
 *                not fixed), 0 elsewhere (double)
 *      C:        0 at fixed points and fixed edges (sinks: sediment reaching
 *                them leaves the model, and their elevation is not changed), 1 elsewhere
 *      Kr, Ks:   bedrock and sediment erodibility, scalars or K x J
 *      Vs:       effective settling velocity, scalar or K x J
 *      par:      [dx dy dt m n phi Ff Hstar tol maxIter pitLimiter doSediment relax adapt]
 *                relax: relaxation factor w for the inflow used for deposition,
 *                QsDep <- w*QsIn + (1-w)*QsDep after each iteration (1 = none).
 *                adapt: if 1, w is multiplied by 0.7 (down to 0.1) whenever the
 *                iteration oscillates: successive changes of Z point in opposite
 *                directions overall (negative dot product) and the largest change
 *                fell by less than half. relax and adapt are optional (default 1, 1).
 *      bvec:     boundary codes [left right upper lower] (0 fixed, 1 mirror, 2 periodic)
 *      QsLat:    optional lateral sediment input to each cell (m^3/yr, solid volume),
 *                added to the cell's inflow from upstream: hillslope soil delivered
 *                to channel cells that do not diffuse (NLDiff), so that deposition
 *                decides how much of it settles and how much travels on. May be
 *                omitted or [] (no lateral input).
 *
 *    Return arguments:
 *      Z, R, H:  values at the end of the step
 *      Qs:       sediment flux leaving each cell (m^3/yr, solid volume)
 *      QsIn:     sediment flux entering each cell
 *      info:     [iterations nLimited produced washload storage exported
 *                 deficit residual converged maxChange relax lateral], with
 *                 residual = produced + lateral - storage - exported + deficit;
 *                 the budget terms are
 *                 solid volumes over the step, weighted by the fraction of
 *                 each cell inside the domain (1/2 on mirror edges, 1/4 at
 *                 mirror corners)
 *      lakeFrac: fraction of each lake's volume filled this step, one value
 *                per lake in order of increasing label
 *      diagnostics (optional): changeHist, the largest change of Z in each
 *                iteration (1 x iterations), and dZlast, the change of Z in
 *                the last iteration (K x J)
 *
 *    Without sediment, R, H, Qs, QsIn and lakeFrac are returned empty and
 *    only the first two entries of info are set.
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2026 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>
#include "flowstack.h"
#include "receivers.h"


/* Solve Z = Zt0 - alpha*S^n for one cell, where the slope is S = lambda*Z - sigma
   (linear in Z once the receivers' new elevations are known). In terms of S:
   S + lambda*alpha*S^n = S0, with S0 = lambda*Zt0 - sigma. The left side
   increases with S, so there is one root in [0, S0], found by Newton's method
   safeguarded by bisection (n = 1 has a closed form). With the pit limiter,
   the cell may not end up below its lowest receiver zmin; *lim is set to 1
   if that happens. */
static double SolveZ(double Zt0, double sigma, double lambda, double alpha, double n,
                     double zmin, int pitLimiter, int *lim) {
    double S0, beta, S, lo, hi, h, dh, Snew, Z;
    int it;

    *lim = 0;
    S0 = lambda*Zt0 - sigma;
    if (S0 <= 0 || lambda == 0) return Zt0; /* no downhill slope: no erosion */
    beta = lambda*alpha;
    if (n == 1) {
        S = S0/(1 + beta);
    } else {
        lo = 0; hi = S0; S = S0;
        for (it = 0; it < 100; it++) {
            h = S + beta*pow(S, n) - S0;
            if (h > 0) hi = S; else lo = S;
            dh = 1 + beta*n*pow(S, n - 1);
            Snew = S - h/dh;
            if (!(Snew > lo && Snew < hi)) Snew = 0.5*(lo + hi);
            if (fabs(Snew - S) <= 1e-14*S0) { S = Snew; break; }
            S = Snew;
        }
    }
    Z = (S + sigma)/lambda;
    if (pitLimiter && Z < zmin) {
        Z = zmin;
        *lim = 1;
    }
    return Z;
}


/* Erosion rate with a smoothed threshold:
   F(x; Ec) = x - Ec*(1 - exp(-x/Ec)), x = K*Omega, Ec the threshold erosion
   rate; F = x for Ec = 0. A series for x/Ec < 1e-3, as in SmoothThreshold.m. */
static double SmoothF(double x, double Ec) {
    double y;
    if (Ec <= 0) return x;
    y = x/Ec;
    if (y < 1e-3) return Ec*y*y*(0.5 - y/6 + y*y/24);
    return x - Ec*(1 - exp(-y));
}

/* dF/dx = 1 - exp(-x/Ec) (1 for Ec = 0) */
static double SmoothdF(double x, double Ec) {
    if (Ec <= 0) return 1;
    return -expm1(-x/Ec);
}

/* F and dF/dx together, with one exponential */
static void SmoothFdF(double x, double Ec, double *F, double *dF) {
    double y, E;
    if (Ec <= 0) { *F = x; *dF = 1; return; }
    y = x/Ec;
    if (y < 1e-3) {
        *F = Ec*y*y*(0.5 - y/6 + y*y/24);
        *dF = y*(1 - y/2 + y*y/6);
        return;
    }
    E = exp(-y);
    *F = x - Ec*(1 - E);
    *dF = 1 - E;
}

/* Solve Z = Zt0 - dt*G(S) for one cell with smoothed thresholds, where
     G(S) = ar*F(Kr*Qm*S^n; Ecr) + as*F(Ks*Qm*S^n; Ecs),  S = lambda*Z - sigma
   (section 14; for Ecr = Ecs = 0 this is SolveZ with alpha = dt*(ar*Kr + as*Ks)*Qm).
   In terms of S: h(S) = S + lambda*dt*G(S) - S0 = 0. If G(S0) > 0 a root lies in
   [0, S0] (h(0) < 0 < h(S0)) and is found by Newton's method safeguarded by
   bisection; it is unique when G increases (ar >= 0). Otherwise (bulking, net
   raising) Z = Zt0 - dt*G with the slope of the current iterate Zcur, as in the
   unthresholded solver. Pit limiter as in SolveZ. */
static long nSolveZT = 0, nNewtonZT = 0; /* diagnostics: calls and Newton iterations (info(13)) */

static double SolveZT(double Zt0, double sigma, double lambda, double dt, double Qm, double n,
                      double ar, double as, double Kr, double Ks, double Ecr, double Ecs,
                      double Zcur, double zmin, int pitLimiter, int *lim) {
    double S0, S, Sn, G, lo, hi, h, dh, Snew, Z, xr, xs;
    int it;

    *lim = 0;
    S0 = lambda*Zt0 - sigma;
    if (S0 <= 0 || lambda == 0) return Zt0; /* no downhill slope: no erosion */
    Sn = pow(S0, n);
    G = ar*SmoothF(Kr*Qm*Sn, Ecr) + as*SmoothF(Ks*Qm*Sn, Ecs);
    if (!(G > 0)) {
        S = lambda*Zcur - sigma;
        Sn = pow(S > 0 ? S : 0, n);
        return Zt0 - dt*(ar*SmoothF(Kr*Qm*Sn, Ecr) + as*SmoothF(Ks*Qm*Sn, Ecs));
    }
    nSolveZT++;
    lo = 0; hi = S0; S = S0;
    /* bounds for n = 1 and ar >= 0: x - Ec <= F(x) <= x gives the root without
       thresholds as a lower bound and the root of the linear form with F = x - Ec
       as an upper bound; h is convex (F is), so Newton started from the upper
       bound converges monotonically */
    if (n == 1 && ar >= 0) {
        double c1 = lambda*dt*(ar*Kr + as*Ks)*Qm;
        double Sup = (S0 + lambda*dt*(ar*Ecr + as*Ecs))/(1 + c1);
        lo = S0/(1 + c1);
        if (Sup < hi) hi = Sup;
        S = hi;
    }
    for (it = 0; it < 100; it++) {
        double Fr = 0, dFr = 0, Fs = 0, dFs = 0, Snm1;
        if (n == 1) { Sn = S; Snm1 = 1; }
        else { Snm1 = pow(S, n - 1); Sn = (S > 0) ? S*Snm1 : 0; }
        xr = Kr*Qm*Sn; xs = Ks*Qm*Sn;
        if (ar != 0) SmoothFdF(xr, Ecr, &Fr, &dFr); /* (as = 0 on bare rock: skip the sediment term) */
        if (as != 0) SmoothFdF(xs, Ecs, &Fs, &dFs);
        h = S + lambda*dt*(ar*Fr + as*Fs) - S0;
        if (fabs(h) <= 1e-12*S0) break; /* residual at round-off level */
        if (h > 0) hi = S; else lo = S;
        dh = 1 + lambda*dt*(ar*Kr*dFr + as*Ks*dFs)*Qm*n*Snm1;
        nNewtonZT++;
        Snew = S - h/dh;
        if (!(Snew > lo && Snew < hi)) Snew = 0.5*(lo + hi);
        if (fabs(Snew - S) <= 1e-12*S0) { S = Snew; break; }
        S = Snew;
    }
    Z = (S + sigma)/lambda;
    if (pitLimiter && Z < zmin) {
        Z = zmin;
        *lim = 1;
    }
    return Z;
}


/* Solve the sediment thickness equation for one cell, with Omega fixed:
     psi(H) = H - H0 - dt/(1-phi)*D(e) + dt*Ks*Om*(1-e) = 0,  e = exp(-H/Hstar),
     D(e) = d0 + f*Om*((1-phi)*Ks*(1-e) + (1-Ff)*Kr*e).
   psi increases with H, psi(0) <= 0 and psi(Hmax) >= 0, so the root is found by
   Newton's method safeguarded by bisection, starting from the estimate Hg. */
static double SolveH(double H0, double Hg, double dt, double Om, double d0, double f, double Ks, double Kr,
                     double phi, double Ff, double Hstar) {
    double Dmax, lo, hi, H, e, D, psi, dpsi, Hnew;
    int it;

    Dmax = d0 + f*Om*(((1 - phi)*Ks > (1 - Ff)*Kr) ? (1 - phi)*Ks : (1 - Ff)*Kr);
    lo = 0; hi = H0 + dt*Dmax/(1 - phi);
    if (hi <= 0) return 0;
    H = Hg < lo ? lo : (Hg > hi ? hi : Hg);
    for (it = 0; it < 100; it++) {
        e = exp(-H/Hstar);
        D = d0 + f*Om*((1 - phi)*Ks*(1 - e) + (1 - Ff)*Kr*e);
        psi = H - H0 - dt*D/(1 - phi) + dt*Ks*Om*(1 - e);
        if (psi > 0) hi = H; else lo = H;
        dpsi = 1 + dt*Om*e/Hstar*(Ks*(1 - f) + f*(1 - Ff)*Kr/(1 - phi));
        Hnew = H - psi/dpsi;
        if (!(Hnew > lo && Hnew < hi)) Hnew = 0.5*(lo + hi);
        if (fabs(Hnew - H) <= 1e-12*(1 + hi)) { H = Hnew; break; }
        H = Hnew;
    }
    return H;
}


/*
 * Shore cells that flow from lake L back into it. Lake surfaces filled by
 * Priority-Flood+epsilon are not quite flat: the filled cells rise by tiny
 * increments away from the spill point. A cell that was not flooded and lies
 * at the lake level can then receive flow from a higher cell of the lake and
 * pass it to a lower one, so the lake drains into a cell upstream of itself.
 * Such cells are added to the lake (lakeOf[c] = L): sediment passes through
 * them like any lake cell and they do not erode in this step. Flow goes
 * strictly downhill on the filled surface Zf, so a path back into the lake
 * stays below the lake's highest cell; the search is confined to cells
 * outside any lake with Zf at most that high. Returns the number of cells
 * added.
 */
static int AbsorbShore(int N, int L, const int *rptr, const int *ridx, const double *Zf,
                       int *lakeOf, int *list, char *mark) {
    double zmax = -HUGE_VAL;
    int c, e, r, n = 0, head = 0, changed, added = 0;

    for (c = 0; c < N; c++) if (lakeOf[c] == L && Zf[c] > zmax) zmax = Zf[c];

    /* cells reached from the lake through cells outside lakes, no higher than zmax */
    for (c = 0; c < N; c++) {
        if (lakeOf[c] != L) continue;
        for (e = rptr[c]; e < rptr[c+1]; e++) {
            r = ridx[e];
            if (lakeOf[r] < 0 && !mark[r] && Zf[r] <= zmax) { mark[r] = 1; list[n++] = r; }
        }
    }
    while (head < n) {
        c = list[head++];
        for (e = rptr[c]; e < rptr[c+1]; e++) {
            r = ridx[e];
            if (lakeOf[r] < 0 && !mark[r] && Zf[r] <= zmax) { mark[r] = 1; list[n++] = r; }
        }
    }

    /* of those, the cells that drain back into the lake (mark = 2) */
    do {
        changed = 0;
        for (head = 0; head < n; head++) {
            c = list[head];
            if (mark[c] == 2) continue;
            for (e = rptr[c]; e < rptr[c+1]; e++) {
                r = ridx[e];
                if (lakeOf[r] == L || mark[r] == 2) { mark[c] = 2; changed = 1; break; }
            }
        }
    } while (changed);

    for (head = 0; head < n; head++) {
        c = list[head];
        if (mark[c] == 2) { lakeOf[c] = L; added++; }
        mark[c] = 0;
    }
    return added;
}


/*
 * Order in which the flux sweep visits cells: every cell after all of its
 * donors (upstream first), and every cell of a lake before any cell outside
 * the lake that receives flow from it. Built by a topological sort (Kahn's
 * algorithm, as for the stack in flowstack.h) of the flow graph with one
 * extra node per lake (numbered N..N+nL-1) between the lake's cells and the
 * cells they drain into. Edges leaving a lake are marked in isExit and listed
 * per lake in exitPtr / exitEdge.
 *
 * If a lake drains into a cell that is upstream of the lake, no such order
 * exists. Shore cells that flow out of the lake and back into it are then
 * added to the lake (AbsorbShore; lakeOf is changed) and the order is built
 * again. If that does not remove the cycle, the lake's extra node is dropped
 * (lakeNode[L] = 0), its sediment passes through untrapped, and a warning is
 * issued. Returns the number of entries in order (N plus the number of lake
 * nodes kept).
 */
static int FluxOrder(int N, int nL, const int *rptr, const int *ridx, const double *Zf, int *lakeOf,
                     char *lakeNode, char *isExit, int *exitPtr, int *exitEdge, int *order) {
    int *indeg, *pos, *list;
    char *done, *mark;
    int M = N + nL, attempt, c, e, L, v, head, tail, nbad, cnt, added;

    indeg = (int *)mxMalloc(M*sizeof(int));
    pos   = (int *)mxMalloc((nL + 1)*sizeof(int));
    done  = (char *)mxMalloc(M);
    list  = NULL; mark = NULL;

    for (attempt = 0; attempt < 6; attempt++) {

        /* mark the edges leaving each lake and list them per lake */
        memset(exitPtr, 0, (nL + 1)*sizeof(int));
        for (c = 0; c < N; c++) {
            L = lakeOf[c];
            for (e = rptr[c]; e < rptr[c+1]; e++) {
                isExit[e] = (L >= 0 && lakeNode[L] && lakeOf[ridx[e]] != L);
                if (isExit[e]) exitPtr[L+1]++;
            }
        }
        for (L = 0; L < nL; L++) exitPtr[L+1] += exitPtr[L];
        memcpy(pos, exitPtr, nL*sizeof(int));
        for (c = 0; c < N; c++)
            for (e = rptr[c]; e < rptr[c+1]; e++)
                if (isExit[e]) exitEdge[pos[lakeOf[c]]++] = e;

        /* in-degrees: an exit edge c -> r becomes c -> lake node -> r */
        memset(indeg, 0, M*sizeof(int));
        for (c = 0; c < N; c++)
            for (e = rptr[c]; e < rptr[c+1]; e++) {
                if (isExit[e]) { indeg[N + lakeOf[c]]++; indeg[ridx[e]]++; }
                else indeg[ridx[e]]++;
            }

        /* Kahn's algorithm */
        head = 0; tail = 0;
        memset(done, 0, M);
        for (v = 0; v < M; v++) if (indeg[v] == 0) order[tail++] = v;
        while (head < tail) {
            v = order[head++];
            done[v] = 1;
            if (v < N) {
                for (e = rptr[v]; e < rptr[v+1]; e++) {
                    int w = isExit[e] ? N + lakeOf[v] : ridx[e];
                    if (--indeg[w] == 0) order[tail++] = w;
                }
            } else {
                L = v - N;
                for (cnt = exitPtr[L]; cnt < exitPtr[L+1]; cnt++) {
                    int w = ridx[exitEdge[cnt]];
                    if (--indeg[w] == 0) order[tail++] = w;
                }
            }
        }

        if (tail == M) {
            /* drop the extra nodes of lakes that have none */
            cnt = 0;
            for (v = 0; v < M; v++)
                if (order[v] < N || lakeNode[order[v] - N]) order[cnt++] = order[v];
            mxFree(indeg); mxFree(pos); mxFree(done);
            if (list) { mxFree(list); mxFree(mark); }
            return cnt;
        }

        /* a cycle: first add to each lake that could not be ordered the shore
           cells that flow out of it and back in, and try again */
        if (!list) {
            list = (int *)mxMalloc(N*sizeof(int));
            mark = (char *)mxCalloc(N, 1);
        }
        added = 0;
        for (L = 0; L < nL; L++)
            if (lakeNode[L] && !done[N + L]) added += AbsorbShore(N, L, rptr, ridx, Zf, lakeOf, list, mark);
        if (added > 0) continue;

        /* otherwise drop the extra nodes of those lakes */
        nbad = 0;
        for (L = 0; L < nL; L++)
            if (lakeNode[L] && !done[N + L]) { lakeNode[L] = 0; nbad++; }
        if (nbad == 0) {
            for (c = 0; c < N; c++)
                if (!done[c] && lakeOf[c] >= 0 && lakeNode[lakeOf[c]]) { lakeNode[lakeOf[c]] = 0; nbad++; }
        }
        mexWarnMsgIdAndTxt("Tadpole:lakeOrder",
            "%d lake(s) drain into cells upstream of themselves; their sediment is not trapped this step.", nbad);
    }
    mxFree(indeg); mxFree(pos); mxFree(done);
    if (list) { mxFree(list); mxFree(mark); }
    mexErrMsgIdAndTxt("Tadpole:lakeOrder", "Could not order the cells for the sediment flux sweep.");
    return 0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    const double *Z, *R = NULL, *H = NULL, *QsPrev = NULL, *W, *StackD, *Q, *Zf, *Lake, *erodes, *Cfix;
    const double *Kr, *Ks = NULL, *Vs = NULL, *par, *bvec, *Lat = NULL;
    double *Zn, *Rn = NULL, *Hn = NULL, *QsOut = NULL, *QsIn = NULL, *info, *lakeFrac = NULL;
    double dx, dy, dt, m, n, phi = 0, Ff = 0, Hstar = 1, tol, a, relax = 1, Ecr = 0, Ecs = 0;
    int adapt = 1;
    int K, J, N, i, j, k, c, t, e, nnz, maxIter, pitLimiter, doSediment, sKr, sKs = 0, sVs = 0;
    int *stack, *rptr, *ridx;
    double *rwsl, *rwacc, *lambda, *Qm;
    int mirror[4];
    long nLimited = 0;

    if (nrhs != 16 && nrhs != 17)
        mexErrMsgIdAndTxt("mexSPACE:nrhs",
            "Usage: [Z, R, H, Qs, QsIn, info, lakeFrac, changeHist, dZlast] = mexSPACE(Z, R, H, QsPrev, W, Stack, Q, Zf, Lake, erodes, C, Kr, Ks, Vs, par, bvec, QsLat)");
    if (nlhs > 9)
        mexErrMsgIdAndTxt("mexSPACE:nlhs", "mexSPACE returns at most 9 outputs.");
    for (i = 0; i < nrhs; i++)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
            mexErrMsgIdAndTxt("mexSPACE:type", "Input %d must be a real double array.", i + 1);
    if (mxGetNumberOfElements(prhs[14]) < 12)
        mexErrMsgIdAndTxt("mexSPACE:par", "par must be [dx dy dt m n phi Ff Hstar tol maxIter pitLimiter doSediment (relax adapt Ecr Ecs)].");

    // Get pointers to inputs
    Z      = mxGetPr(prhs[0]);
    W      = mxGetPr(prhs[4]);
    StackD = mxGetPr(prhs[5]);
    Q      = mxGetPr(prhs[6]);
    Zf     = mxGetPr(prhs[7]);
    Lake   = mxGetPr(prhs[8]);
    erodes = mxGetPr(prhs[9]);
    Cfix   = mxGetPr(prhs[10]);
    Kr     = mxGetPr(prhs[11]);
    par    = mxGetPr(prhs[14]);
    bvec   = mxGetPr(prhs[15]);

    K = (int)mxGetM(prhs[0]);
    J = (int)mxGetN(prhs[0]);
    N = K*J;
    if ((int)mxGetNumberOfElements(prhs[4]) != 8*N || (int)mxGetNumberOfElements(prhs[5]) != N)
        mexErrMsgIdAndTxt("mexSPACE:size", "W must be K x J x 8 and Stack must have K*J elements.");

    dx = par[0]; dy = par[1]; dt = par[2]; m = par[3]; n = par[4];
    tol = par[8]; maxIter = (int)par[9]; pitLimiter = (par[10] != 0); doSediment = (par[11] != 0);
    if (mxGetNumberOfElements(prhs[14]) > 12) relax = par[12];
    if (mxGetNumberOfElements(prhs[14]) > 13) adapt = (par[13] != 0);
    /* smoothed erosion thresholds (m/yr), bedrock and sediment (section 14); 0: none */
    if (mxGetNumberOfElements(prhs[14]) > 14) Ecr = par[14];
    if (mxGetNumberOfElements(prhs[14]) > 15) Ecs = par[15];
    sKr = mxGetNumberOfElements(prhs[11]) > 1; /* stride: 0 for a scalar, 1 for a matrix */

    if (doSediment) {
        R     = mxGetPr(prhs[1]);
        H     = mxGetPr(prhs[2]);
        QsPrev = mxGetPr(prhs[3]);
        Ks    = mxGetPr(prhs[12]);
        Vs    = mxGetPr(prhs[13]);
        sKs = mxGetNumberOfElements(prhs[12]) > 1;
        sVs = mxGetNumberOfElements(prhs[13]) > 1;
        phi = par[5]; Ff = par[6]; Hstar = par[7];
        if ((int)mxGetNumberOfElements(prhs[1]) != N || (int)mxGetNumberOfElements(prhs[2]) != N
            || (int)mxGetNumberOfElements(prhs[3]) != N)
            mexErrMsgIdAndTxt("mexSPACE:size", "R, H and QsPrev must be K x J.");
        if (nrhs > 16 && !mxIsEmpty(prhs[16])) {
            if ((int)mxGetNumberOfElements(prhs[16]) != N)
                mexErrMsgIdAndTxt("mexSPACE:size", "QsLat must be K x J or [].");
            Lat = mxGetPr(prhs[16]);
        }
    }

    // Create arrays for return arguments
    Zn = mxGetPr(plhs[0] = mxCreateDoubleMatrix(K, J, mxREAL));
    memcpy(Zn, Z, N*sizeof(double));
    if (doSediment) {
        Rn    = mxGetPr(plhs[1] = mxCreateDoubleMatrix(K, J, mxREAL));
        Hn    = mxGetPr(plhs[2] = mxCreateDoubleMatrix(K, J, mxREAL));
        QsOut = mxGetPr(plhs[3] = mxCreateDoubleMatrix(K, J, mxREAL));
        QsIn  = mxGetPr(plhs[4] = mxCreateDoubleMatrix(K, J, mxREAL));
        memcpy(Rn, R, N*sizeof(double));
        memcpy(Hn, H, N*sizeof(double));
    } else {
        for (i = 1; i <= 4; i++) plhs[i] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    info = mxGetPr(plhs[5] = mxCreateDoubleMatrix(1, 13, mxREAL));
    nSolveZT = 0; nNewtonZT = 0;

    for (i = 0; i < 4; i++) mirror[i] = (bvec[i] == FS_MIRROR);
    a = dx*dy;

    stack = (int *)mxMalloc(N*sizeof(int));
    for (t = 0; t < N; t++) stack[t] = (int)StackD[t] - 1;

    /*
     * Compact receiver lists (CSR, receivers.h): the receivers of cell c are
     * ridx[rptr[c] .. rptr[c+1]-1], with rwsl the weight / distance for the slope,
     * S = sum rwsl*(Z(c) - Z(receiver)) (mirror-edge weights restored), and rwacc
     * the weight for passing a flux downstream (with the mirror-edge doubling).
     */
    nnz = FSReceiverLists(K, J, W, bvec, dx, dy, &rptr, &ridx, &rwacc, &rwsl, NULL);
    lambda = (double *)mxMalloc(N*sizeof(double)); /* S = lambda*Z(c) - sum rwsl*Z(receivers) */
    Qm     = (double *)mxMalloc(N*sizeof(double));
    for (c = 0; c < N; c++) {
        lambda[c] = 0;
        for (e = rptr[c]; e < rptr[c+1]; e++) lambda[c] += rwsl[e];
        Qm[c] = (m == 1) ? Q[c] : pow(Q[c], m);
    }

    if (!doSediment) {

        /* ---------- stream power only: one sweep, receivers first ----------
           A river flowing into a lake erodes down to the lake's water surface
           (the filled surface), not to the lake floor, because stream power
           depends on the water surface slope. Lake cells do not erode, so their
           surface stays fixed. (With flood = 0 there are no lakes and Zf = Z.) */
        for (c = 0; c < N; c++) if (Lake[c] > 0) Zn[c] = Zf[c];
        for (t = 0; t < N; t++) {
            double sigma = 0, zmin = HUGE_VAL, zr;
            int lim;
            c = stack[t];
            if (erodes[c] == 0 || rptr[c+1] == rptr[c]) continue;
            for (e = rptr[c]; e < rptr[c+1]; e++) {
                zr = Zn[ridx[e]];
                sigma += rwsl[e]*zr;
                if (zr < zmin) zmin = zr;
            }
            if (Ecr > 0)
                Zn[c] = SolveZT(Z[c], sigma, lambda[c], dt, Qm[c], n, 1, 0, Kr[sKr*c], 0, Ecr, 0, Z[c], zmin, pitLimiter, &lim);
            else
                Zn[c] = SolveZ(Z[c], sigma, lambda[c], dt*Kr[sKr*c]*Qm[c], n, zmin, pitLimiter, &lim);
            nLimited += lim;
        }
        for (c = 0; c < N; c++) if (Lake[c] > 0) Zn[c] = Z[c];
        info[0] = 1;
        info[1] = (double)nLimited;
        for (i = 6; i < nlhs; i++) plhs[i] = mxCreateDoubleMatrix(0, 0, mxREAL);

    } else {

        /* ---------- stream power with sediment ---------- */
        double *fdep, *d0c, *om, *Zold, *QsNew, *QsDep, *dZ, *lakeVol, *trapped, *exitQ;
        int *lakeOf, *labelLake, *order, *exitPtr, *exitEdge, nL = 0, nOrder, iter, converged = 0, L;
        char *lakeNode, *isExit;
        double exported = 0, deficit = 0, maxChange = 0, prevChange = 0, dot;
        double *hist = (nlhs > 7) ? (double *)mxMalloc((maxIter + 1)*sizeof(double)) : NULL;
        double produced = 0, washload = 0, storage = 0, lateral = 0;

        /* coefficients of the local flux formula: Qs_out = (Qs_in + a*prod)/(1 + V*a/Q), D = V*Qs_out/Q */
        fdep = (double *)mxMalloc(N*sizeof(double)); /* fraction of local production deposited in the cell */
        d0c  = (double *)mxMalloc(N*sizeof(double)); /* D = d0c*Qs_in + fdep*prod */
        for (c = 0; c < N; c++) {
            double V = Vs[sVs*c];
            if (Q[c] > 0) {
                fdep[c] = (V*a/Q[c])/(1 + V*a/Q[c]);
                d0c[c]  = (V/Q[c])/(1 + V*a/Q[c]);
            } else {
                fdep[c] = 1; d0c[c] = 1/a; /* no water: everything arriving is deposited */
            }
        }

        /* fraction of each cell inside the reflected domain */
        om = (double *)mxMalloc(N*sizeof(double));
        for (j = 0; j < J; j++)
            for (i = 0; i < K; i++) {
                double f = 1;
                if (mirror[0] && j == 0)   f /= 2;
                if (mirror[1] && j == J-1) f /= 2;
                if (mirror[2] && i == 0)   f /= 2;
                if (mirror[3] && i == K-1) f /= 2;
                om[i + K*j] = f;
            }

        /* lakes: each lake is labelled by its spill cell; number them 0..nL-1 in
           order of increasing label. Lake volume is weighted by om so that it is
           consistent with the fluxes. */
        labelLake = (int *)mxMalloc((N + 1)*sizeof(int));
        for (c = 0; c <= N; c++) labelLake[c] = -1;
        for (c = 0; c < N; c++) if (Lake[c] > 0) labelLake[(int)Lake[c]] = 0;
        for (c = 1; c <= N; c++) if (labelLake[c] == 0) labelLake[c] = nL++;
        lakeOf = (int *)mxMalloc(N*sizeof(int));
        for (c = 0; c < N; c++) lakeOf[c] = (Lake[c] > 0) ? labelLake[(int)Lake[c]] : -1;
        lakeVol  = (double *)mxCalloc(nL + 1, sizeof(double));
        trapped  = (double *)mxCalloc(nL + 1, sizeof(double));
        lakeNode = (char *)mxMalloc(nL + 1);
        for (c = 0; c < N; c++)
            if (lakeOf[c] >= 0 && Zf[c] > Z[c]) lakeVol[lakeOf[c]] += om[c]*a*(Zf[c] - Z[c]);
        for (L = 0; L < nL; L++) lakeNode[L] = 1;
        plhs[6] = mxCreateDoubleMatrix(nL, 1, mxREAL);
        lakeFrac = mxGetPr(plhs[6]);

        /* order for the flux sweep: upstream first, with every cell of a lake
           before any cell that the lake drains into, so that a lake's total
           sediment inflow is known when the trapping is computed */
        order    = (int *)mxMalloc((N + nL)*sizeof(int));
        isExit   = (char *)mxMalloc(nnz + 1);
        exitPtr  = (int *)mxMalloc((nL + 1)*sizeof(int));
        exitEdge = (int *)mxMalloc((nnz + 1)*sizeof(int));
        exitQ    = (double *)mxCalloc(nnz + 1, sizeof(double)); /* flux along each edge leaving a lake, before trapping */
        if (nL == 0) {
            for (t = 0; t < N; t++) order[t] = stack[N - 1 - t];
            nOrder = N;
            memset(isExit, 0, nnz + 1);
            exitPtr[0] = 0;
        } else {
            nOrder = FluxOrder(N, nL, rptr, ridx, Zf, lakeOf, lakeNode, isExit, exitPtr, exitEdge, order);
        }

        Zold  = (double *)mxMalloc(N*sizeof(double));
        QsNew = (double *)mxMalloc(N*sizeof(double)); /* flux entering each cell */
        dZ    = (double *)mxCalloc(N, sizeof(double));  /* change of Z in the last iteration */

        /* inflow used for deposition (under-relaxed between iterations). Starting
           guess: the previous step's outflow routed through the current weights, so
           that cells whose donors changed with the routing start from their new
           donors' fluxes, plus the lateral input */
        QsDep = (double *)mxCalloc(N, sizeof(double));
        for (c = 0; c < N; c++)
            for (e = rptr[c]; e < rptr[c+1]; e++)
                QsDep[ridx[e]] += rwacc[e]*QsPrev[c];
        if (Lat) for (c = 0; c < N; c++) QsDep[c] += Lat[c];

        for (iter = 1; iter <= maxIter; iter++) {

            memcpy(Zold, Zn, N*sizeof(double));

            /* ---------- (a) receivers-first sweep: new elevations ----------
               A river flowing into a lake erodes down to the lake's water surface
               (its base level), not to the lake floor, so receivers in a lake count
               with the elevation of the filled surface, which is fixed during the
               step. Other receivers have already been updated in this sweep. */
            for (t = 0; t < N; t++) {
                double d0, sigma, zmin, zr, Zt0, ee, Znew = 0, Hnew = 0, Rnew, Om, Hg;
                double Krc, Ksc;
                int it, lim = 0;
                c = stack[t];
                if (Cfix[c] == 0 || lakeOf[c] >= 0) continue; /* lake cells are updated after the flux sweep */
                d0 = d0c[c]*QsDep[c]; /* deposition from incoming sediment */
                if (erodes[c] == 0 || rptr[c+1] == rptr[c]) {
                    /* deposition only */
                    double dz = dt*d0/(1 - phi);
                    Zn[c] = Z[c] + dz; Hn[c] = H[c] + dz; Rn[c] = R[c];
                    continue;
                }
                Krc = Kr[sKr*c]; Ksc = Ks[sKs*c];
                sigma = 0; zmin = HUGE_VAL;
                for (e = rptr[c]; e < rptr[c+1]; e++) {
                    int r = ridx[e];
                    zr = (lakeOf[r] >= 0) ? Zf[r] : Zn[r];
                    sigma += rwsl[e]*zr;
                    if (zr < zmin) zmin = zr;
                }
                Zt0 = Z[c] + dt*d0/(1 - phi);
                if (Ecr > 0 || Ecs > 0) {
                    /* no stream power even at the steepest slope the cell can have in this
                       step (S0, at Zt0): deposition only (and both cell solves are skipped) */
                    double S0c = lambda[c]*Zt0 - sigma, w0 = Qm[c]*pow(S0c > 0 ? S0c : 0, n);
                    if (Krc*w0 == 0 && Ksc*w0 == 0) {
                        double dz = dt*d0/(1 - phi);
                        Zn[c] = Z[c] + dz; Hn[c] = H[c] + dz; Rn[c] = R[c];
                        continue;
                    }
                }
                /* The cover factor e and H depend on each other. A few passes per cell
                   are enough, because e is refreshed in every Gauss-Seidel iteration,
                   so at convergence the full system for Z and H is satisfied. */
                ee = exp(-Hn[c]/Hstar);
                Hg = Hn[c];
                for (it = 0; it < 3; it++) {
                    double Kc   = Krc*ee + Ksc*(1 - ee);                     /* bed lowering per unit Omega */
                    double kap  = (1 - phi)*Ksc*(1 - ee) + (1 - Ff)*Krc*ee;  /* solid production per unit Omega */
                    double ceff = Kc - fdep[c]*kap/(1 - phi);                /* net lowering per unit Omega */
                    double enew, S;
                    if (Ecr > 0 || Ecs > 0) {
                        /* smoothed thresholds (section 14): Er = e F(Kr Om; Ecr) and
                           Es = (1 - e) F(Ks Om; Ecs); the H equation takes the rates
                           F(Ks Om) and F(Kr Om) in place of Ks Om and Kr Om */
                        double ar = ee*(1 - fdep[c]*(1 - Ff)/(1 - phi));
                        double as = (1 - ee)*(1 - fdep[c]);
                        Znew = SolveZT(Zt0, sigma, lambda[c], dt, Qm[c], n, ar, as, Krc, Ksc, Ecr, Ecs,
                                       Zn[c], zmin, pitLimiter, &lim);
                        S = lambda[c]*Znew - sigma;
                        Om = Qm[c]*pow(S > 0 ? S : 0, n);
                        Hnew = SolveH(H[c], Hg, dt, 1.0, d0, fdep[c], SmoothF(Ksc*Om, Ecs), SmoothF(Krc*Om, Ecr), phi, Ff, Hstar);
                    } else {
                        if (ceff > 0) {
                            Znew = SolveZ(Zt0, sigma, lambda[c], dt*ceff*Qm[c], n, zmin, pitLimiter, &lim);
                        } else {
                            /* locally redeposited rock-derived sediment bulks up more than the
                               bed is lowered: use Omega from the current iterate */
                            S = lambda[c]*Zn[c] - sigma;
                            Znew = Zt0 - dt*ceff*Qm[c]*pow(S > 0 ? S : 0, n);
                            lim = 0;
                        }
                        S = lambda[c]*Znew - sigma;
                        Om = Qm[c]*pow(S > 0 ? S : 0, n);
                        Hnew = SolveH(H[c], Hg, dt, Om, d0, fdep[c], Ksc, Krc, phi, Ff, Hstar);
                    }
                    enew = exp(-Hnew/Hstar);
                    Hg = Hnew;
                    if (fabs(enew - ee) < 1e-12) break;
                    ee = enew;
                }
                Rnew = Znew - Hnew;
                if (lim) {
                    /* the pit limiter raised the cell: keep the bedrock no higher than before */
                    if (Rnew > R[c]) { Rnew = R[c]; Hnew = Znew - Rnew; }
                    nLimited++;
                }
                Zn[c] = Znew; Hn[c] = Hnew; Rn[c] = Rnew;
            }

            /* ---------- (b) upstream-first sweep: sediment fluxes ----------
               Each cell's inflow starts with its lateral input (soil delivered by
               hillslopes), to which the fluxes from its donors are added */
            if (Lat) memcpy(QsNew, Lat, N*sizeof(double));
            else memset(QsNew, 0, N*sizeof(double));
            exported = 0; deficit = 0;
            for (t = 0; t < nOrder; t++) {
                int v = order[t];
                double qin, qout;
                if (v >= N) {
                    /* all cells of lake L have been visited: trap a fraction tau of the
                       sediment leaving it (up to the lake's volume), pass on the rest */
                    double Qin = 0, tau = 0;
                    int cnt;
                    L = v - N;
                    for (cnt = exitPtr[L]; cnt < exitPtr[L+1]; cnt++) {
                        e = exitEdge[cnt];
                        Qin += om[ridx[e]]*exitQ[e]; /* om-weighted total inflow to the lake */
                    }
                    if (Qin > 0) {
                        tau = (1 - phi)*lakeVol[L]/(Qin*dt);
                        if (tau > 1) tau = 1;
                    }
                    trapped[L] = tau*Qin*dt;
                    for (cnt = exitPtr[L]; cnt < exitPtr[L+1]; cnt++) {
                        e = exitEdge[cnt];
                        QsNew[ridx[e]] += (1 - tau)*exitQ[e];
                    }
                    continue;
                }
                c = v;
                qin = QsNew[c];
                QsIn[c] = qin;
                if (Cfix[c] == 0) {
                    /* fixed point or fixed edge: sediment leaves the model */
                    exported += qin*om[c]*dt;
                    continue;
                }
                if (lakeOf[c] >= 0) {
                    qout = qin; /* sediment passes through lake cells; trapping happens at the outlet */
                } else {
                    qout = qin + a*(-(1 - phi)*(Hn[c] - H[c]) - (1 - Ff)*(Rn[c] - R[c]))/dt;
                    if (qout < 0) {
                        deficit -= qout*om[c]*dt;
                        qout = 0;
                    }
                }
                QsOut[c] = qout;
                if (rptr[c+1] == rptr[c]) {
                    /* an outlet that is not a fixed point (with no fixed boundaries or
                       sinks, the routing drains the grid through its lowest cell): the
                       sediment leaving it leaves the model */
                    exported += qout*om[c]*dt;
                }
                for (e = rptr[c]; e < rptr[c+1]; e++) {
                    double q = rwacc[e]*qout;
                    if (isExit[e]) exitQ[e] = q;
                    else QsNew[ridx[e]] += q;
                }
            }

            /* raise every lake cell by the same fraction of its water depth */
            for (L = 0; L < nL; L++)
                lakeFrac[L] = (lakeVol[L] > 0) ? trapped[L]/((1 - phi)*lakeVol[L]) : 0;
            for (c = 0; c < N; c++) {
                if (lakeOf[c] >= 0) {
                    double dz = lakeFrac[lakeOf[c]]*((Zf[c] > Z[c]) ? Zf[c] - Z[c] : 0);
                    Zn[c] = Z[c] + dz; Hn[c] = H[c] + dz; Rn[c] = R[c];
                }
            }

            /* ---------- (c) convergence ---------- */
            maxChange = 0; dot = 0;
            for (c = 0; c < N; c++) {
                double d = Zn[c] - Zold[c];
                dot += d*dZ[c];
                dZ[c] = d;
                if (fabs(d) > maxChange) maxChange = fabs(d);
            }
            if (hist) hist[iter-1] = maxChange;
            if (maxChange < tol) { converged = 1; break; }

            /* under-relaxation: if the iteration oscillates (the change of Z
               reverses overall and is not decaying fast), reduce w */
            if (adapt && iter >= 3 && dot < 0 && maxChange > 0.5*prevChange) {
                relax *= 0.7;
                if (relax < 0.1) relax = 0.1;
            }
            prevChange = maxChange;
            for (c = 0; c < N; c++) QsDep[c] = relax*QsIn[c] + (1 - relax)*QsDep[c];
        }
        if (iter > maxIter) iter = maxIter;
        if (nlhs > 7) {
            double *h = mxGetPr(plhs[7] = mxCreateDoubleMatrix(1, iter, mxREAL));
            memcpy(h, hist, iter*sizeof(double));
        }
        if (nlhs > 8) {
            double *d = mxGetPr(plhs[8] = mxCreateDoubleMatrix(K, J, mxREAL));
            memcpy(d, dZ, N*sizeof(double));
        }
        if (hist) mxFree(hist);

        /* mass budget of this step (solid volumes, om-weighted) */
        for (c = 0; c < N; c++) {
            produced += om[c]*a*(1 - Ff)*(R[c] - Rn[c]); /* sediment made from bedrock by rivers */
            washload += om[c]*a*Ff*(R[c] - Rn[c]);       /* bedrock lost as wash load */
            storage  += om[c]*a*(1 - phi)*(Hn[c] - H[c]); /* change in stored sediment */
            if (Lat) lateral += om[c]*Lat[c]*dt;          /* lateral input (soil from hillslopes) */
        }
        info[0] = iter;
        info[1] = (double)nLimited;
        info[2] = produced;
        info[3] = washload;
        info[4] = storage;
        info[5] = exported;
        info[6] = deficit;
        info[7] = produced + lateral - storage - exported + deficit;
        info[8] = converged;
        info[9] = maxChange;
        info[10] = relax;
        info[11] = lateral;
        info[12] = nSolveZT ? (double)nNewtonZT/nSolveZT : 0; /* mean Newton iterations per thresholded cell solve */

        mxFree(fdep); mxFree(d0c); mxFree(om); mxFree(labelLake); mxFree(lakeOf);
        mxFree(lakeVol); mxFree(trapped); mxFree(lakeNode); mxFree(order); mxFree(isExit);
        mxFree(exitPtr); mxFree(exitEdge); mxFree(exitQ); mxFree(Zold); mxFree(QsNew);
        mxFree(QsDep); mxFree(dZ);
    }

    mxFree(stack); mxFree(rptr); mxFree(ridx); mxFree(rwsl); mxFree(rwacc);
    mxFree(lambda); mxFree(Qm);
}
