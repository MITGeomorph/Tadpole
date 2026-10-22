/*
 * =============================================================
 * flowstack.h - drainage accumulation and stack ordering, shared by
 *               mexD8, mexDinf and mexDms
 *
 * FlowAccumulate() takes a K x J x 8 array of flow routing weights W
 * and computes drainage area and upstream-weighted sums, without
 * recursion. It also returns, if requested, the stack (an ordering of
 * all cells in which every cell comes before all of its donors) and
 * the receiver of each nonzero weight.
 *
 * The cells are ordered with a topological sort (Kahn's algorithm),
 * as in the multiple-flow-direction stack of FastScape: count each
 * cell's donors, start from the cells that have none, and add a cell
 * to the order once all of its donors have been added. This visits
 * every cell before any of its receivers, for single (D8) or multiple
 * (D-infinity, multi-slope) flow directions, in O(N) time. The stack
 * is that order reversed, so the cells at its top are the outlets, as
 * in Braun & Willett (2013).
 *
 * Weight directions k = 0..7 point to the neighbor at
 *   (i + di[k], j + dj[k]), di = {0,-1,-1,-1,0,1,1,1}, dj = {1,1,0,-1,-1,-1,0,1}
 * wrapped across periodic boundaries. Weights that point off the grid
 * across a fixed or mirror boundary are ignored.
 *
 * Mirror boundaries: a cell on a mirror boundary receiving flow from a
 * cell off that boundary receives twice the donor's weighted area,
 * because the donor's mirror image also drains to it.
 *
 * This is a header for MATLAB MEX-files written in C or C++.
 * =============================================================
 */

#ifndef FLOWSTACK_H
#define FLOWSTACK_H

#include "mex.h"

#define FS_MIRROR   1
#define FS_PERIODIC 2

static const int FS_DI[8] = {0, -1, -1, -1, 0, 1, 1, 1};
static const int FS_DJ[8] = {1, 1, 0, -1, -1, -1, 0, 1};


/* Linear index of the neighbor of (i,j) in direction k, or -1 if it is off the grid */
static int FSNeighbor(int i, int j, int k, int K, int J, const double bdy[]) {
    int ni = i + FS_DI[k], nj = j + FS_DJ[k];
    if (ni < 0)       { if (bdy[2] != FS_PERIODIC) return -1; ni += K; }
    else if (ni >= K) { if (bdy[3] != FS_PERIODIC) return -1; ni -= K; }
    if (nj < 0)       { if (bdy[0] != FS_PERIODIC) return -1; nj += J; }
    else if (nj >= J) { if (bdy[1] != FS_PERIODIC) return -1; nj -= J; }
    return ni + K * nj;
}


/* Factor applied to flow that receiver (i,j) gets from its neighbor in direction k:
   doubled for each mirror boundary the receiver is on and the donor is not (the
   donor's mirror image across that boundary also drains to the receiver). A corner
   between two mirror boundaries receiving from a diagonal neighbor gets a factor 4. */
static double FSMirrorFactor(int i, int j, int k, int K, int J, const double bdy[]) {
    double f = 1;
    if (i == 0 && bdy[2] == FS_MIRROR && (k == 5 || k == 6 || k == 7)) f *= 2;
    else if (i == K-1 && bdy[3] == FS_MIRROR && (k == 1 || k == 2 || k == 3)) f *= 2;
    if (j == 0 && bdy[0] == FS_MIRROR && (k == 7 || k == 0 || k == 1)) f *= 2;
    else if (j == J-1 && bdy[1] == FS_MIRROR && (k == 3 || k == 4 || k == 5)) f *= 2;
    return f;
}


/*
 * K, J:     number of rows and columns
 * W:        K x J x 8 flow routing weights (column-major)
 * bdy:      boundary codes [left right upper lower] (0=fixed, 1=mirror, 2=periodic)
 * P:        K x J values to sum upstream (precipitation rates)
 * A:        output, K x J drainage area in cells
 * Pav:      output, K x J upstream weighted sum of P. Divide by A to get
 *           upstream averages.
 * Stack:    output or NULL, K*J 1-based linear indices, each cell before all of its donors
 * Receiver: output or NULL, K x J x 8 1-based linear indices of the receiver of
 *           each nonzero weight (0 elsewhere); must be zeroed by the caller
 */
static void FlowAccumulate(int K, int J, const double W[], const double bdy[],
                           const double P[], double A[], double Pav[],
                           double Stack[], double Receiver[]) {
    const int N = K * J;
    int *ndon, *order, *rcv;
    int c, i, j, k, r, head = 0, tail = 0;

    ndon  = (int *)mxCalloc(N, sizeof(int));     /* number of donors not yet visited */
    order = (int *)mxMalloc(N * sizeof(int));    /* cells in upstream-to-downstream order */
    rcv   = (int *)mxMalloc(8 * (size_t)N * sizeof(int));  /* receiver in each direction, or -1 */

    /* find receivers and count donors */
    for (j = 0; j < J; j++) {
        for (i = 0; i < K; i++) {
            c = i + K * j;
            for (k = 0; k < 8; k++) {
                r = (W[k*N + c] > 0) ? FSNeighbor(i, j, k, K, J, bdy) : -1;
                rcv[8*c + k] = r;
                if (r >= 0) {
                    ndon[r]++;
                    if (Receiver) Receiver[k*N + c] = r + 1;
                }
            }
            A[c] = 1;  /* each cell drains at least itself */
            Pav[c] = P[c];
        }
    }

    /* start from the cells that have no donors */
    for (c = 0; c < N; c++) {
        if (ndon[c] == 0) order[tail++] = c;
    }

    /* pass flow downstream; a cell joins the order when all of its donors have */
    while (head < tail) {
        c = order[head++];
        for (k = 0; k < 8; k++) {
            r = rcv[8*c + k];
            if (r < 0) continue;
            {
                const int ir = r % K, jr = r / K;
                const double w = W[k*N + c] * FSMirrorFactor(ir, jr, (k + 4) % 8, K, J, bdy);
                A[r] += w * A[c];
                Pav[r] += w * Pav[c];
            }
            if (--ndon[r] == 0) order[tail++] = r;
        }
    }

    if (tail != N) {
        mxFree(ndon); mxFree(order); mxFree(rcv);
        mexErrMsgIdAndTxt("flowstack:loop", "The flow routing weights contain a loop, so drainage area is undefined.");
    }

    if (Stack) {
        for (c = 0; c < N; c++) Stack[N - 1 - c] = order[c] + 1;
    }

    mxFree(ndon);
    mxFree(order);
    mxFree(rcv);
}

#endif /* FLOWSTACK_H */
