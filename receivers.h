/*
 * =============================================================
 * receivers.h - compact receiver lists built from the flow routing
 *               weights, used by mexSPACE
 *
 * The routing functions (mexD8, mexDinf, mexDms) return K x J x 8
 * arrays of weights W. Solvers that sweep the stack need, for each
 * cell, only its receivers, so FSReceiverLists() stores them in
 * compressed sparse row (CSR) form: the receivers of cell c (0-based
 * linear index) are ridx[rptr[c] .. rptr[c+1]-1], in direction order
 * k = 0..7, with cells in column-major order.
 *
 * For each receiver it stores:
 *   rwacc: weight for passing a flux downstream, including the doubling
 *          of flow that a mirror-boundary cell receives from the interior
 *          (FSMirrorFactor), as in the drainage area calculation.
 *   rwsl:  (optional) weight / distance for the slope,
 *          S = sum rwsl*(Z(c) - Z(receiver)). On mirror boundaries the
 *          routing functions halve the weights of flow from a boundary cell
 *          into the interior (FixMirrorWeights); those weights are doubled
 *          again here, as in UpwindSlope and ReceiverLists.
 *   rinvL: (optional) 1/distance to the receiver.
 *
 * Pass NULL for an optional output that is not needed. The arrays are
 * allocated with mxMalloc (with one spare element, so they are never of
 * size zero); the caller frees them. Returns the number of receivers.
 *
 * This is a header for MATLAB MEX-files written in C.
 * =============================================================
 */

#ifndef RECEIVERS_H
#define RECEIVERS_H

#include <math.h>
#include "mex.h"
#include "flowstack.h"

static int FSReceiverLists(int K, int J, const double W[], const double bdy[], double dx, double dy,
                           int **rptr, int **ridx, double **rwacc, double **rwsl, double **rinvL) {
    const int N = K * J;
    const double diagd = sqrt(dx*dx + dy*dy);
    double invd[8];
    int mirror[4];
    int i, j, k, c, e, nnz = 0;

    invd[0] = 1/dx; invd[1] = 1/diagd; invd[2] = 1/dy; invd[3] = 1/diagd;
    invd[4] = 1/dx; invd[5] = 1/diagd; invd[6] = 1/dy; invd[7] = 1/diagd;
    for (i = 0; i < 4; i++) mirror[i] = (bdy[i] == FS_MIRROR);

    for (c = 0; c < N; c++)
        for (k = 0; k < 8; k++)
            if (W[k*N + c] > 0) nnz++;
    *rptr  = (int *)mxMalloc((N + 1)*sizeof(int));
    *ridx  = (int *)mxMalloc((nnz + 1)*sizeof(int));
    *rwacc = (double *)mxMalloc((nnz + 1)*sizeof(double));
    if (rwsl)  *rwsl  = (double *)mxMalloc((nnz + 1)*sizeof(double));
    if (rinvL) *rinvL = (double *)mxMalloc((nnz + 1)*sizeof(double));

    e = 0;
    for (j = 0; j < J; j++) {
        for (i = 0; i < K; i++) {
            c = i + K*j;
            (*rptr)[c] = e;
            for (k = 0; k < 8; k++) {
                double w = W[k*N + c], fs = 1;
                int r;
                if (w <= 0) continue;
                r = FSNeighbor(i, j, k, K, J, bdy);
                if (r < 0) continue;
                (*ridx)[e]  = r;
                (*rwacc)[e] = w*FSMirrorFactor(r % K, r / K, (k + 4) % 8, K, J, bdy);
                if (rwsl) {
                    /* directions whose weights FixMirrorWeights halved on each mirror boundary */
                    if (mirror[0] && j == 0   && (k == 7 || k == 0 || k == 1)) fs *= 2;
                    if (mirror[1] && j == J-1 && (k == 3 || k == 4 || k == 5)) fs *= 2;
                    if (mirror[2] && i == 0   && (k == 5 || k == 6 || k == 7)) fs *= 2;
                    if (mirror[3] && i == K-1 && (k == 1 || k == 2 || k == 3)) fs *= 2;
                    (*rwsl)[e] = invd[k]*(fs*w);
                }
                if (rinvL) (*rinvL)[e] = invd[k];
                e++;
            }
        }
    }
    (*rptr)[N] = e;
    return e;
}

#endif /* RECEIVERS_H */
