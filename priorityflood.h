/*
 * =============================================================
 * priorityflood.h - depression filling by Priority-Flood+epsilon,
 *                   shared by mexFlood, mexD8, mexDinf and mexDms
 *
 * PriorityFloodEpsilon() fills depressions and flats in an elevation
 * grid so that every cell has at least one D8 neighbor that is
 * strictly lower, except sinks and cells on fixed boundaries. Any
 * flow routing scheme that only sends flow to strictly lower
 * neighbors (D8, D-infinity, multi-slope) therefore produces flow
 * paths that end at sinks or fixed boundaries and contain no loops.
 *
 * Algorithm: Priority-Flood+epsilon (Barnes, Lehman & Mulla, 2014,
 * Computers & Geosciences 62:117-127, Algorithm 3), with D8
 * connectivity and periodic/mirror boundaries. O(N log N), and O(N)
 * when the input has no local minima (the grid is then returned
 * unchanged, because Priority-Flood+epsilon would not modify it).
 *
 * Boundary codes (bvec = [left right upper lower]):
 *   0 = fixed:    cells on that edge are sinks.
 *   1 = mirror:   the grid is reflected across that edge. A reflected
 *                 neighbor is always also an ordinary neighbor of the
 *                 same cell, so mirror edges add no connections.
 *   2 = periodic: the grid wraps across that edge. Periodic edges
 *                 must be paired (left with right, upper with lower).
 *
 * Elevation increment: a filled cell is raised to at least
 *   max(nextafter(z), z + inc)
 * above the cell it drains to, where inc depends on epsilon:
 *   epsilon <  0 (PF_EPS_AUTO): inc = 1e-12 * max|Z|, but at least
 *                 1e-100, so that elevation differences in filled areas
 *                 never approach the floating-point underflow limit
 *                 (which breaks D-infinity, whose slopes are squared,
 *                 and multi-slope, which divides by the sum of slopes).
 *   epsilon == 0: inc = 0, i.e. the smallest representable increment
 *                 (std::nextafter), as in RichDEM.
 *   epsilon >  0: inc = epsilon.
 *
 * Work buffers are kept between calls (in each MEX file that includes
 * this header) and are freed when the MEX file is cleared.
 *
 * This is a header for MATLAB MEX-files written in C or C++.
 * =============================================================
 */

#ifndef PRIORITYFLOOD_H
#define PRIORITYFLOOD_H

#include "mex.h"
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#define PF_FIXED    0
#define PF_MIRROR   1
#define PF_PERIODIC 2

#define PF_EPS_AUTO (-1.0)
#define PF_EPS_REL  1.0e-12
#define PF_EPS_MIN  1.0e-100

/* bits in the per-cell flag array */
#define PF_CLOSED 1   /* cell has received its final elevation */
#define PF_EDGE   2   /* cell is on the edge of the grid */

typedef struct {
    double z;
    int idx;
} PFCell;

typedef struct {
    int K, J;                 /* rows, columns */
    int bl, br, bu, bd;       /* boundary codes: left, right, upper, lower */
    int off[8];               /* linear index offsets for interior cells */
} PFGrid;

/* D8 neighbor offsets (row, column) */
static const int PF_DI[8] = {-1, 1,  0, 0, -1, -1,  1, 1};
static const int PF_DJ[8] = { 0, 0, -1, 1, -1,  1, -1, 1};

/* Work buffers, kept between calls (each cell enters each of them at most once) */
static unsigned char *pf_flags = NULL;
static int           *pf_pit   = NULL;
static PFCell        *pf_heap  = NULL;
static size_t         pf_cap   = 0;


/* ===== WORK BUFFERS ===== */

static void PFFreeBuffers(void) {
    if (pf_flags) mxFree(pf_flags);
    if (pf_pit)   mxFree(pf_pit);
    if (pf_heap)  mxFree(pf_heap);
    pf_flags = NULL;
    pf_pit = NULL;
    pf_heap = NULL;
    pf_cap = 0;
}

static void PFReserve(size_t N) {
    unsigned char *f;
    int *p;
    PFCell *h;

    if (N <= pf_cap) return;
    PFFreeBuffers();
    /* if any allocation fails, MATLAB frees the others and the statics stay NULL */
    f = (unsigned char *)mxMalloc(N * sizeof(unsigned char));
    p = (int *)mxMalloc(N * sizeof(int));
    h = (PFCell *)mxMalloc(N * sizeof(PFCell));
    mexMakeMemoryPersistent(f);
    mexMakeMemoryPersistent(p);
    mexMakeMemoryPersistent(h);
    pf_flags = f;
    pf_pit = p;
    pf_heap = h;
    pf_cap = N;
    mexAtExit(PFFreeBuffers);
}


/* ===== MIN-HEAP ON ELEVATION ===== */

static void PFHeapPush(PFCell *h, size_t *n, double z, int idx) {
    size_t pos = (*n)++;
    while (pos > 0) {
        size_t parent = (pos - 1) / 2;
        if (h[parent].z <= z) break;
        h[pos] = h[parent];
        pos = parent;
    }
    h[pos].z = z;
    h[pos].idx = idx;
}

static int PFHeapPop(PFCell *h, size_t *n) {
    int top = h[0].idx;
    size_t sz = --(*n), pos = 0;
    PFCell last = h[sz];
    if (sz == 0) return top;
    while (1) {
        size_t child = 2 * pos + 1;
        if (child >= sz) break;
        if (child + 1 < sz && h[child + 1].z < h[child].z) child++;
        if (h[child].z >= last.z) break;
        h[pos] = h[child];
        pos = child;
    }
    h[pos] = last;
    return top;
}


/* ===== NEIGHBORS ===== */

/* Neighbors of cell c; returns their number. Neighbors across fixed or mirror
   edges are omitted (a mirrored neighbor is a duplicate of another neighbor
   of the same cell). */
static int PFNeighbors(const PFGrid *g, int c, int edge, int nb[8]) {
    int k, n = 0, i, j, ni, nj;
    if (!edge) {
        for (k = 0; k < 8; k++) nb[k] = c + g->off[k];
        return 8;
    }
    i = c % g->K;
    j = c / g->K;
    for (k = 0; k < 8; k++) {
        ni = i + PF_DI[k];
        nj = j + PF_DJ[k];
        if (ni < 0)          { if (g->bu != PF_PERIODIC) continue; ni += g->K; }
        else if (ni >= g->K) { if (g->bd != PF_PERIODIC) continue; ni -= g->K; }
        if (nj < 0)          { if (g->bl != PF_PERIODIC) continue; nj += g->J; }
        else if (nj >= g->J) { if (g->br != PF_PERIODIC) continue; nj -= g->J; }
        nb[n++] = ni + g->K * nj;
    }
    return n;
}


/* ===== PRIORITY-FLOOD+EPSILON ===== */

/*
 * K, J:    number of rows and columns
 * Z:       K*J input elevations (column-major); not modified
 * Zf:      K*J output: filled elevations (must not overlap Z)
 * sinksD, sinksL: K*J sink mask as double or as logical (nonzero = sink);
 *          pass at most one of them, or neither (NULL) for no sinks
 * bvec:    boundary codes [left right upper lower]
 * epsilon: elevation increment; see the top of this file
 * minima:  K*J output, or NULL. If given, must be zeroed by the caller, and
 *          is set to 1 where a non-sink cell of the INPUT surface has no
 *          strictly lower neighbor.
 *
 * lake: output or NULL, K x J, zeroed by the caller. Each cell raised by
 *       filling gets the 1-based linear index of the spill cell of its
 *       depression (the unraised cell it ultimately drains to), so all cells
 *       of one lake share a label. Cells that are not raised stay 0.
 *
 * Returns 1 if the input contained local minima (i.e. Zf may differ from Z),
 * 0 otherwise. If there are no sinks and no fixed boundaries, the lowest cell
 * is used as the only outlet and a warning (mexFlood:noSinks) is issued.
 */
static int PriorityFloodEpsilon(int K, int J, const double *Z, double *Zf,
                                const double *sinksD, const mxLogical *sinksL,
                                const double *bvec, double epsilon, mxLogical *minima, double *lake) {
    PFGrid g;
    int b[4], k, i, j, c, m, n, nb[8], anyMinima = 0;
    size_t N, nHeap = 0, pitHead = 0, pitTail = 0, cc;
    unsigned char *flg;
    double inc, maxAbs = 0.0;

    /* ----- check arguments ----- */
    if (K <= 0 || J <= 0)
        mexErrMsgIdAndTxt("mexFlood:Z", "The elevation grid must not be empty.");
    if ((double)K * (double)J >= (double)INT_MAX)
        mexErrMsgIdAndTxt("mexFlood:Z", "The elevation grid has too many elements.");
    N = (size_t)K * (size_t)J;

    for (k = 0; k < 4; k++) {
        if (bvec[k] != PF_FIXED && bvec[k] != PF_MIRROR && bvec[k] != PF_PERIODIC)
            mexErrMsgIdAndTxt("mexFlood:bvec", "Boundary codes must be 0 (fixed), 1 (mirror) or 2 (periodic).");
        b[k] = (int)bvec[k];
    }
    if ((b[0] == PF_PERIODIC) != (b[1] == PF_PERIODIC) || (b[2] == PF_PERIODIC) != (b[3] == PF_PERIODIC))
        mexErrMsgIdAndTxt("mexFlood:bvec", "Periodic boundaries must be paired (left with right, upper with lower).");
    if ((b[0] == PF_MIRROR || b[1] == PF_MIRROR) && J < 2)
        mexErrMsgIdAndTxt("mexFlood:bvec", "Mirror left/right boundaries require at least 2 columns.");
    if ((b[2] == PF_MIRROR || b[3] == PF_MIRROR) && K < 2)
        mexErrMsgIdAndTxt("mexFlood:bvec", "Mirror upper/lower boundaries require at least 2 rows.");

    for (cc = 0; cc < N; cc++) {
        double a = fabs(Z[cc]);
        if (!(a <= DBL_MAX))   /* NaN or Inf */
            mexErrMsgIdAndTxt("mexFlood:Z", "Elevations must be finite (no NaN or Inf).");
        if (a > maxAbs) maxAbs = a;
    }
    if (epsilon < 0.0) {
        inc = PF_EPS_REL * maxAbs;
        if (inc < PF_EPS_MIN) inc = PF_EPS_MIN;
    } else {
        inc = epsilon;
    }

    g.K = K; g.J = J;
    g.bl = b[0]; g.br = b[1]; g.bu = b[2]; g.bd = b[3];
    for (k = 0; k < 8; k++) g.off[k] = PF_DI[k] + K * PF_DJ[k];

    PFReserve(N);
    flg = pf_flags;

    /* ----- mark edges, sinks and fixed-boundary cells; seed the priority queue ----- */
    for (j = 0; j < J; j++) {
        for (i = 0; i < K; i++) {
            int edge, seed;
            c = i + K * j;
            edge = (i == 0 || i == K - 1 || j == 0 || j == J - 1);
            seed = edge && ((j == 0 && g.bl == PF_FIXED) || (j == J - 1 && g.br == PF_FIXED) ||
                            (i == 0 && g.bu == PF_FIXED) || (i == K - 1 && g.bd == PF_FIXED));
            if (!seed) seed = sinksL ? (sinksL[c] != 0) : (sinksD && sinksD[c] != 0.0);
            flg[c] = (unsigned char)((edge ? PF_EDGE : 0) | (seed ? PF_CLOSED : 0));
            if (seed) PFHeapPush(pf_heap, &nHeap, Z[c], c);
        }
    }

    /* ----- local minima of the input surface (PF_CLOSED marks exactly the sinks here) ----- */
    for (c = 0; c < (int)N; c++) {
        double zc;
        int lower = 0;
        if (flg[c] & PF_CLOSED) continue;
        n = PFNeighbors(&g, c, flg[c] & PF_EDGE, nb);
        zc = Z[c];
        for (k = 0; k < n; k++) {
            if (Z[nb[k]] < zc) { lower = 1; break; }
        }
        if (!lower) {
            anyMinima = 1;
            if (!minima) break;
            minima[c] = 1;
        }
    }

    memcpy(Zf, Z, N * sizeof(double));

    /* With no minima, every steepest-descent path already ends at a sink, and
       Priority-Flood+epsilon would leave every cell unchanged. */
    if (!anyMinima && nHeap > 0) return 0;

    if (nHeap == 0) {
        /* nowhere for flow to go: use the lowest cell as the single outlet */
        int cmin = 0;
        for (c = 1; c < (int)N; c++) if (Z[c] < Z[cmin]) cmin = c;
        flg[cmin] |= PF_CLOSED;
        PFHeapPush(pf_heap, &nHeap, Z[cmin], cmin);
        mexWarnMsgIdAndTxt("mexFlood:noSinks",
            "No fixed boundaries or sinks were given; the lowest cell is used as the outlet.");
    }

    /* ----- flood ----- */
    while (1) {
        double zc, zraise;
        if (pitHead < pitTail) {
            c = pf_pit[pitHead++];
        } else if (nHeap > 0) {
            c = PFHeapPop(pf_heap, &nHeap);
        } else {
            break;
        }

        zc = Zf[c];
        zraise = nextafter(zc, HUGE_VAL);
        if (zc + inc > zraise) zraise = zc + inc;

        n = PFNeighbors(&g, c, flg[c] & PF_EDGE, nb);
        for (k = 0; k < n; k++) {
            m = nb[k];
            if (flg[m] & PF_CLOSED) continue;
            flg[m] |= PF_CLOSED;
            if (Zf[m] <= zc) {
                /* m is in a depression or flat: raise it so it drains to c */
                Zf[m] = zraise;
                if (lake) lake[m] = (lake[c] > 0.0) ? lake[c] : (double)(c + 1);
                pf_pit[pitTail++] = m;
            } else {
                PFHeapPush(pf_heap, &nHeap, Zf[m], m);
            }
        }
    }

    return anyMinima;
}

#endif /* PRIORITYFLOOD_H */
