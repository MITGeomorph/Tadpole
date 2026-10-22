/*
 * =============================================================
 * hilldiffusion.h - implicit diffusion steps on the conservative 5-point
 *                   operator of AndrewsHanks / HillslopeFaces, shared by
 *                   mexHillLinear.c and mexHillPicard.c
 *                  
 *
 *    The operator, for a cell diffusivity D (zero at cells that don't diffuse):
 *
 *      RHS_c = sum over the four faces f of c_f*(z_nbr - z_c),
 *      c_f = (D_c + D_nbr)/(2*h_f^2)
 *
 *    On a face between a diffusing cell and a cell that doesn't diffuse (a
 *    sink), the diffusing cell's D is used instead of the average, so the flux
 *    into sinks is not halved. Neighbours across mirror edges are the reflected
 *    cells, across periodic edges the wrapped cells, and beyond fixed edges
 *    ghosts of zero elevation and zero D (as in PadGrid). Sinks keep their
 *    elevation.
 *
 *    dodiff: 1 at diffusing cells, 0 at sinks, and 2 at channel cells that do
 *    not diffuse: sinks that keep their elevation
 *    and receive material from higher diffusing neighbours, but give none to
 *    lower ones: a face between a channel cell and a diffusing neighbour is
 *    closed (no flux) where the channel cell is the higher. The comparison
 *    uses the reference elevations Zref of HDSolveBE (Picard passes its
 *    current iterate; mexHillLinear repeats the solve until the faces agree
 *    with the solution). Only HDSolveBE supports dodiff = 2.
 *
 *    HDFaces       face neighbours, face conductances, and the fraction of
 *                  each cell inside the domain (omega)
 *    HDSolveADI    one Peaceman-Rachford ADI step
 *    HDSolveBE     one backward Euler step by IC(0)-preconditioned conjugate
 *                  gradients (workspace HDWork, allocated once by HDAlloc)
 *    Both compute the new elevations and the delivery to sinks from the face
 *    fluxes, so the volume budget closes to round-off.
 *
 * Copyright (c) 2026 Taylor Perron
 * =============================================================
 */

#include <math.h>
#include <string.h>


/* Neighbour of cell (i,j) across face f (0 right, 1 left, 2 up, 3 down), as in
   PadGrid: the reflected cell across a mirror edge, the wrapped cell across a
   periodic edge, -1 (a ghost) beyond a fixed edge. */
static int FaceNeighbor(int i, int j, int f, int K, int J, const double *bvec) {
    switch (f) {
        case 0:
            if (j < J-1) return i + K*(j+1);
            if (bvec[1] == 1) return i + K*(J-2);
            if (bvec[1] == 2) return i;
            return -1;
        case 1:
            if (j > 0) return i + K*(j-1);
            if (bvec[0] == 1) return i + K;
            if (bvec[0] == 2) return i + K*(J-1);
            return -1;
        case 2:
            if (i > 0) return i-1 + K*j;
            if (bvec[2] == 1) return 1 + K*j;
            if (bvec[2] == 2) return K-1 + K*j;
            return -1;
        default:
            if (i < K-1) return i+1 + K*j;
            if (bvec[3] == 1) return K-2 + K*j;
            if (bvec[3] == 2) return K*j;
            return -1;
    }
}


/* Face neighbours nbr (4*N), face conductances cf (4*N) for the cell
   diffusivity D (D[sD*c]; sD = 0 for a scalar), and omega (N). */
static void HDFaces(int K, int J, const double *dodiff, const double *D, int sD, double dx, double dy,
                    const double *bvec, int *nbr, double *cf, double *om) {
    double h2[4];
    int i, j, f, c;
    h2[0] = dx*dx; h2[1] = dx*dx; h2[2] = dy*dy; h2[3] = dy*dy;
    for (j = 0; j < J; j++) {
        for (i = 0; i < K; i++) {
            double w = 1;
            c = i + K*j;
            if (bvec[0] == 1 && j == 0)   w /= 2;
            if (bvec[1] == 1 && j == J-1) w /= 2;
            if (bvec[2] == 1 && i == 0)   w /= 2;
            if (bvec[3] == 1 && i == K-1) w /= 2;
            om[c] = w;
            for (f = 0; f < 4; f++) {
                int n = FaceNeighbor(i, j, f, K, J, bvec);
                int diffc = dodiff[c] == 1, diffn = 0, sinkn = 0;
                double Dc, Dn = 0, sum;
                Dc = diffc ? D[sD*c] : 0;
                if (n >= 0) {
                    diffn = dodiff[n] == 1;
                    sinkn = !diffn;
                    Dn = diffn ? D[sD*n] : 0;
                }
                sum = Dc + Dn;
                /* one-sided face between a diffusing cell and a sink: the diffusing cell's D */
                if ((diffc && sinkn) || (!diffc && diffn)) cf[4*c+f] = sum/h2[f];
                else cf[4*c+f] = 0.5*sum/h2[f];
                nbr[4*c+f] = n;
            }
        }
    }
}


/* One step of the Thomas algorithm (forward elimination) for element k of a
   line of length len, at linear index idx (prev = index of element k-1). The
   line's matrix is I - h*L, with L given by its sub-diagonal a, diagonal d and
   super-diagonal b (couplings to elements k-1 and k+1, wrapped on periodic
   lines). x holds the right side on input. With cyclic lines (Sherman-Morrison,
   as in Numerical Recipes' cyclic), the corners are removed from the matrix, the
   first and last diagonal entries modified, and a second system with right
   side u = (gam, 0, ..., 0, alp) is solved in zz. */
static void ElimStep(int k, int len, int idx, int prev, int cyclic, double gam, double alp, double bet,
                     const double *a, const double *d, const double *b, double h,
                     double *cp, double *x, double *zz) {
    double sub = -h*a[idx], dia = 1 - h*d[idx], den;
    if (cyclic && k == 0) dia -= gam;
    if (cyclic && k == len-1) dia -= alp*bet/gam;
    den = (k > 0) ? dia - sub*cp[prev] : dia;
    cp[idx] = (k < len-1) ? -h*b[idx]/den : 0;
    x[idx] = ((k > 0) ? x[idx] - sub*x[prev] : x[idx])/den;
    if (cyclic) {
        double u = (k == 0) ? gam : ((k == len-1) ? alp : 0);
        zz[idx] = ((k > 0) ? u - sub*zz[prev] : u)/den;
    }
}


/* Solve (I - h*L) x = r along nl lines of length len (element k of line m at
   index m*sl + k*ss); r is in x on input. Lines along x (ss = K) are swept all
   together, element by element, so memory is accessed contiguously; lines
   along y (ss = 1) one at a time. Periodic lines are cyclic tridiagonal. */
static void SolveLines(int nl, int len, int sl, int ss, int cyclic, const double *a, const double *d,
                       const double *b, double h, double *x, double *cp, double *zz, double *fact) {
    int m, k, idx, first, last;
    double gam, alp, bet;
#define LINE_ENDS(m) first = (m)*sl; last = (m)*sl + (len-1)*ss; \
    gam = -(1 - h*d[first]); alp = -h*b[last]; bet = -h*a[first];
    if (ss == 1) {
        for (m = 0; m < nl; m++) {
            LINE_ENDS(m)
            for (k = 0; k < len; k++) {
                idx = m*sl + k;
                ElimStep(k, len, idx, idx - 1, cyclic, gam, alp, bet, a, d, b, h, cp, x, zz);
            }
            for (k = len-2; k >= 0; k--) {
                idx = m*sl + k;
                x[idx] -= cp[idx]*x[idx + 1];
                if (cyclic) zz[idx] -= cp[idx]*zz[idx + 1];
            }
            if (cyclic) {
                double fc = (x[first] + bet*x[last]/gam)/(1 + zz[first] + bet*zz[last]/gam);
                for (k = 0; k < len; k++) x[m*sl + k] -= fc*zz[m*sl + k];
            }
        }
    } else {
        for (k = 0; k < len; k++)
            for (m = 0; m < nl; m++) {
                LINE_ENDS(m)
                idx = m*sl + k*ss;
                ElimStep(k, len, idx, idx - ss, cyclic, gam, alp, bet, a, d, b, h, cp, x, zz);
            }
        for (k = len-2; k >= 0; k--)
            for (m = 0; m < nl; m++) {
                idx = m*sl + k*ss;
                x[idx] -= cp[idx]*x[idx + ss];
                if (cyclic) zz[idx] -= cp[idx]*zz[idx + ss];
            }
        if (cyclic) {
            for (m = 0; m < nl; m++) {
                LINE_ENDS(m)
                fact[m] = (x[first] + bet*x[last]/gam)/(1 + zz[first] + bet*zz[last]/gam);
            }
            for (k = 0; k < len; k++)
                for (m = 0; m < nl; m++) x[m*sl + k*ss] -= fact[m]*zz[m*sl + k*ss];
        }
    }
#undef LINE_ENDS
}


/* One Peaceman-Rachford ADI step (Crank-Nicolson in 2D), with RHS split into
   its x and y faces, Lx + Ly:
     (I - dt/2 Lx) u* = (I + dt/2 Ly) z^n,   (I - dt/2 Ly) u = (I + dt/2 Lx) u*
   Lx, Ly by lines: for each diffusing cell, a = coupling to the previous
   element of its line, b = to the next (wrapped on periodic lines; a reflected
   neighbour across a mirror edge is an interior element of the line), d =
   -(sum of the conductances of its two faces). Sinks have zero rows, so their
   rows of I - h*L are identity rows. New elevations from the face fluxes:
   z^n + dt*Lx u* + dt/2*Ly (z^n + u), which equals u to round-off; sinks
   receive the same expression (Deliv). */
static void HDSolveADI(int K, int J, const double *Z, const double *dodiff, const int *nbr, const double *cf,
                       const double *bvec, double dt, double *Zn, double *Deliv) {
    int N = K*J, i, j, c, f, pos;
    double *xa, *xd, *xb, *ya, *yd, *yb, *us, *u, *cp, *zz, *fact, hh = 0.5*dt;
    for (c = 0; c < N; c++)
        if (dodiff[c] == 2)
            mexErrMsgIdAndTxt("hilldiffusion:dodiff", "The ADI step does not support channel cells that lose material (dodiff = 2).");
    xa = (double *)mxCalloc(N, sizeof(double)); xd = (double *)mxCalloc(N, sizeof(double));
    xb = (double *)mxCalloc(N, sizeof(double)); ya = (double *)mxCalloc(N, sizeof(double));
    yd = (double *)mxCalloc(N, sizeof(double)); yb = (double *)mxCalloc(N, sizeof(double));
    for (j = 0; j < J; j++)
        for (i = 0; i < K; i++) {
            c = i + K*j;
            if (dodiff[c] == 0) continue;
            for (f = 0; f < 4; f++) {
                double v = cf[4*c+f];
                /* position of the neighbour in the line: -1 previous, 1 next, 0 ghost */
                switch (f) {
                    case 0:  pos = (j < J-1) ? 1 : (bvec[1] == 1 ? -1 : (bvec[1] == 2 ? 1 : 0)); break;
                    case 1:  pos = (j > 0) ? -1 : (bvec[0] == 1 ? 1 : (bvec[0] == 2 ? -1 : 0)); break;
                    case 2:  pos = (i > 0) ? -1 : (bvec[2] == 1 ? 1 : (bvec[2] == 2 ? -1 : 0)); break;
                    default: pos = (i < K-1) ? 1 : (bvec[3] == 1 ? -1 : (bvec[3] == 2 ? 1 : 0)); break;
                }
                if (f < 2) {
                    xd[c] -= v;
                    if (pos == -1) xa[c] += v; else if (pos == 1) xb[c] += v;
                } else {
                    yd[c] -= v;
                    if (pos == -1) ya[c] += v; else if (pos == 1) yb[c] += v;
                }
            }
        }
    us   = (double *)mxMalloc(N*sizeof(double));
    u    = (double *)mxMalloc(N*sizeof(double));
    cp   = (double *)mxMalloc(N*sizeof(double));
    zz   = (double *)mxMalloc(N*sizeof(double));
    fact = (double *)mxMalloc((K > J ? K : J)*sizeof(double));

    /* half step 1: (I - h Lx) u* = (I + h Ly) z^n */
    for (j = 0; j < J; j++)
        for (i = 0; i < K; i++) {
            c = i + K*j;
            us[c] = Z[c] + hh*(yd[c]*Z[c] + ya[c]*Z[j*K + (i > 0 ? i-1 : K-1)] + yb[c]*Z[j*K + (i < K-1 ? i+1 : 0)]);
        }
    SolveLines(K, J, 1, K, bvec[0] == 2, xa, xd, xb, hh, us, cp, zz, fact);
    /* half step 2: (I - h Ly) u = (I + h Lx) u* */
    for (j = 0; j < J; j++)
        for (i = 0; i < K; i++) {
            c = i + K*j;
            u[c] = us[c] + hh*(xd[c]*us[c] + xa[c]*us[i + K*(j > 0 ? j-1 : J-1)] + xb[c]*us[i + K*(j < J-1 ? j+1 : 0)]);
        }
    SolveLines(J, K, K, 1, bvec[2] == 2, ya, yd, yb, hh, u, cp, zz, fact);

    for (c = 0; c < N; c++) {
        double fx = 0, fy = 0;
        for (f = 0; f < 4; f++) {
            int n = nbr[4*c+f];
            double v = cf[4*c+f];
            if (v == 0) continue;
            if (f < 2) fx += v*(((n >= 0) ? us[n] : 0) - us[c]);
            else fy += v*(((n >= 0) ? Z[n] + u[n] : 0) - Z[c] - u[c]);
        }
        if (dodiff[c] != 0) {
            Zn[c] = Z[c] + dt*fx + hh*fy;
        } else {
            Zn[c] = Z[c];
            Deliv[c] = dt*fx + hh*fy;
        }
    }
    mxFree(xa); mxFree(xd); mxFree(xb); mxFree(ya); mxFree(yd); mxFree(yb);
    mxFree(us); mxFree(u); mxFree(cp); mxFree(zz); mxFree(fact);
}


/* Workspace for HDSolveBE */
typedef struct {
    int *rptr, *col, *lptr;
    double *val, *diag, *b, *Lval, *Ldiag, *x, *r, *z, *p, *q;
} HDWork;

static void HDAlloc(HDWork *w, int N) {
    w->rptr  = (int *)mxMalloc((N + 1)*sizeof(int));
    w->col   = (int *)mxMalloc(4*N*sizeof(int));
    w->lptr  = (int *)mxMalloc(N*sizeof(int));
    w->val   = (double *)mxMalloc(4*N*sizeof(double));
    w->diag  = (double *)mxMalloc(N*sizeof(double));
    w->b     = (double *)mxMalloc(N*sizeof(double));
    w->Lval  = (double *)mxMalloc((4*N + 1)*sizeof(double));
    w->Ldiag = (double *)mxMalloc(N*sizeof(double));
    w->x = (double *)mxMalloc(N*sizeof(double));
    w->r = (double *)mxMalloc(N*sizeof(double));
    w->z = (double *)mxMalloc(N*sizeof(double));
    w->p = (double *)mxMalloc(N*sizeof(double));
    w->q = (double *)mxMalloc(N*sizeof(double));
}

static void HDFree(HDWork *w) {
    mxFree(w->rptr); mxFree(w->col); mxFree(w->lptr); mxFree(w->val); mxFree(w->diag); mxFree(w->b);
    mxFree(w->Lval); mxFree(w->Ldiag); mxFree(w->x); mxFree(w->r); mxFree(w->z); mxFree(w->p); mxFree(w->q);
}


/* Whether the face between a diffusing cell c and a channel cell n (dodiff = 2)
   is closed: the channel cell, fixed at its start-of-step elevation Z[n], is
   higher than c in the elevations Zr. No hillslope flux crosses a closed face
   */
static int HDClosed(int c, int n, const double *dodiff, const double *Z, const double *Zr) {
    return dodiff[c] == 1 && dodiff[n] == 2 && Z[n] > Zr[c];
}


/* Number of faces between channel cells (dodiff = 2) and diffusing cells whose
   closed state differs between the elevations Za and Zb */
static int HDClosedChanged(int N, const int *nbr, const double *dodiff, const double *Z,
                           const double *Za, const double *Zb) {
    int c, f, n, count = 0;
    for (c = 0; c < N; c++) {
        if (dodiff[c] != 1) continue;
        for (f = 0; f < 4; f++) {
            n = nbr[4*c+f];
            if (n >= 0 && dodiff[n] == 2 && HDClosed(c, n, dodiff, Z, Za) != HDClosed(c, n, dodiff, Z, Zb)) count++;
        }
    }
    return count;
}


/* One backward Euler step, z - dt*RHS(z) = z^n (Z), for the face conductances
   cf. The unknowns are all cells; rows of sinks are identity rows, and their
   coupling to diffusing cells is moved to the right side. Each row is
   multiplied by omega, which makes the matrix symmetric (a mirror-edge cell
   couples to its interior neighbour across two faces): a symmetric M-matrix,
   solved by conjugate gradients with an IC(0) preconditioner on its own
   sparsity pattern (periodic links included), starting from Zg (or Z if NULL).
   The iteration stops when the largest residual per unit cell area is below
   tol times the largest change of z, or below 1e-13 of the largest |z^n|. New
   elevations Zn = z^n + dt*RHS(z) at diffusing cells; sinks receive dt*RHS(z)
   (Deliv). Returns the number of CG iterations; *converged and *maxF set. */
static int HDSolveBE(HDWork *w, int N, const double *Z, const double *Zg, const double *Zref, const double *dodiff,
                     const int *nbr, const double *cf, const double *om, double dt, double tol, int maxIter,
                     double *Zn, double *Deliv, int *converged, double *maxFout) {
    int *rptr = w->rptr, *col = w->col, *lptr = w->lptr;
    double *val = w->val, *diag = w->diag, *b = w->b, *Lval = w->Lval, *Ldiag = w->Ldiag;
    double *x = w->x, *r = w->r, *z = w->z, *p = w->p, *q = w->q;
    double rz, rzold, alpha, beta, pq, maxF, maxDz, zscale;
    int c, f, e, iter;
    /* elevations that decide which faces next to channel cells (dodiff = 2) are
       closed: Zref, else the starting guess Zg, else z^n */
    const double *Zr = Zref ? Zref : (Zg ? Zg : Z);

    *converged = 0;

    /* Matrix rows (omega-scaled): diagonal, and off-diagonal entries to diffusing
       neighbours, merged where two faces reach the same cell and sorted by column.
       Right side: omega*z^n plus the coupling to sinks (ghosts have zero elevation). */
    e = 0;
    for (c = 0; c < N; c++) {
        rptr[c] = e;
        if (dodiff[c] != 1) { /* sinks, including channel cells */
            diag[c] = 1; b[c] = Z[c];
            continue;
        }
        diag[c] = 1; b[c] = Z[c];
        for (f = 0; f < 4; f++) {
            int n = nbr[4*c+f], s, t;
            double a = dt*cf[4*c+f];
            if (n >= 0 && HDClosed(c, n, dodiff, Z, Zr)) continue; /* closed face: no flux */
            diag[c] += a;
            if (n < 0) continue;
            if (dodiff[n] != 1) { b[c] += a*Z[n]; continue; } /* sinks at their elevation */
            for (s = rptr[c]; s < e; s++) if (col[s] == n) break;
            if (s < e) { val[s] -= a; continue; }
            /* insert, keeping the row sorted by column */
            for (t = e; t > rptr[c] && col[t-1] > n; t--) { col[t] = col[t-1]; val[t] = val[t-1]; }
            col[t] = n; val[t] = -a; e++;
        }
        diag[c] *= om[c]; b[c] *= om[c];
        for (f = rptr[c]; f < e; f++) val[f] *= om[c];
    }
    rptr[N] = e;

    /* IC(0) preconditioner, A ~ L*L', L with the sparsity pattern of the lower part
       of A: L(i,j) = (A(i,j) - sum_k L(i,k)*L(j,k))/L(j,j) for j < i,
       L(i,i) = sqrt(A(i,i) - sum_j L(i,j)^2). Rows are sorted, so the lower part
       of row i is its first entries, rptr[i] .. lptr[i]-1. */
    for (c = 0; c < N; c++) {
        int s, s2, t2;
        double d;
        for (s = rptr[c]; s < rptr[c+1] && col[s] < c; s++) ;
        lptr[c] = s;
        d = diag[c];
        for (s = rptr[c]; s < lptr[c]; s++) {
            int jj = col[s];
            double v = val[s];
            /* subtract L(c,k)*L(jj,k) over the common lower columns k < jj */
            for (s2 = rptr[c]; s2 < s; s2++)
                for (t2 = rptr[jj]; t2 < lptr[jj]; t2++)
                    if (col[t2] == col[s2]) { v -= Lval[s2]*Lval[t2]; break; }
            Lval[s] = v/Ldiag[jj];
            d -= Lval[s]*Lval[s];
        }
        Ldiag[c] = (d > 0) ? sqrt(d) : sqrt(diag[c]); /* cannot happen for an M-matrix */
    }

    /* preconditioned conjugate gradients */
    zscale = 0;
    for (c = 0; c < N; c++) {
        x[c] = (Zg && dodiff[c] == 1) ? Zg[c] : Z[c];
        if (fabs(Z[c]) > zscale) zscale = fabs(Z[c]);
    }

    /* initial residual r = b - A*x, and the convergence measures: the largest
       residual per unit area and the largest change of z */
    maxF = 0; maxDz = 0;
    for (c = 0; c < N; c++) {
        int s;
        double Ax = diag[c]*x[c], F, dz;
        for (s = rptr[c]; s < rptr[c+1]; s++) Ax += val[s]*x[col[s]];
        r[c] = b[c] - Ax;
        F = fabs(r[c])/om[c]; dz = fabs(x[c] - Z[c]);
        if (F > maxF) maxF = F;
        if (dz > maxDz) maxDz = dz;
    }

    /* Each iteration makes five passes over the grid (memory traffic dominates
       the cost): forward solve, backward solve with the dot product r'z,
       update of p, A*p with p'Ap, and the update of x and r with the
       convergence measures. */
    rzold = 0;
    for (iter = 0; ; iter++) {
        int s;
        if (maxF <= tol*maxDz || maxF <= 1e-13*zscale) { *converged = 1; break; }
        if (iter == maxIter) break;

        /* z = (L*L')^-1 r: forward solve by rows, backward solve by columns of L */
        for (c = 0; c < N; c++) {
            double v = r[c];
            for (s = rptr[c]; s < lptr[c]; s++) v -= Lval[s]*z[col[s]];
            z[c] = v/Ldiag[c];
        }
        rz = 0;
        for (c = N-1; c >= 0; c--) {
            z[c] /= Ldiag[c]; /* final: all later rows have been subtracted */
            rz += r[c]*z[c];
            for (s = rptr[c]; s < lptr[c]; s++) z[col[s]] -= Lval[s]*z[c];
        }

        if (iter == 0) {
            memcpy(p, z, N*sizeof(double));
        } else {
            beta = rz/rzold;
            for (c = 0; c < N; c++) p[c] = z[c] + beta*p[c];
        }
        rzold = rz;

        pq = 0;
        for (c = 0; c < N; c++) {
            double Ap = diag[c]*p[c];
            for (s = rptr[c]; s < rptr[c+1]; s++) Ap += val[s]*p[col[s]];
            q[c] = Ap;
            pq += p[c]*Ap;
        }
        if (pq <= 0) break; /* only if the residual is already at round-off */
        alpha = rz/pq;
        maxF = 0; maxDz = 0;
        for (c = 0; c < N; c++) {
            double F, dz;
            x[c] += alpha*p[c];
            r[c] -= alpha*q[c];
            F = fabs(r[c])/om[c]; dz = fabs(x[c] - Z[c]);
            if (F > maxF) maxF = F;
            if (dz > maxDz) maxDz = dz;
        }
    }
    if (!*converged)
        mexWarnMsgIdAndTxt("Tadpole:hillLinear",
            "Linear hillslope solver did not converge in %d iterations (residual %g m, change %g m).", iter, maxF, maxDz);

    /* new elevations and delivery to sinks from the face fluxes at the solution;
       closed faces next to channel cells (dodiff = 2) carry no flux */
    for (c = 0; c < N; c++) {
        double rhs = 0;
        for (f = 0; f < 4; f++) {
            int n = nbr[4*c+f];
            double zn = (n >= 0) ? x[n] : 0;
            if (cf[4*c+f] == 0) continue;
            if (n >= 0 && (HDClosed(c, n, dodiff, Z, Zr) || HDClosed(n, c, dodiff, Z, Zr))) continue;
            rhs += cf[4*c+f]*(zn - x[c]);
        }
        if (dodiff[c] == 1) {
            Zn[c] = Z[c] + dt*rhs;
        } else {
            Zn[c] = Z[c];
            Deliv[c] = dt*rhs;
        }
    }
    *maxFout = maxF;
    return iter;
}
