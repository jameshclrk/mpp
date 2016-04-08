/* * MPP Coursework - MPI Edge Reconstruction
 * Copyright (C) 2015,1016 James Clark
 *
 * This file is part of MPP Coursework.
 *
 * MPP Coursework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPP Coursework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPP Coursework.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file serial/update.c
 * @author James Clark
 * @brief Serial Update Code
 */

#include <mpi.h>
#include <math.h>

#include <precision.h>
#include <functions.h>

step_return update_tick (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** edge, real ** old, real ** new){
    int i, j;
    real delta = 0.0;
    step_return retval = {0.0, 0.0};

    /* periodic boundary conditions for left and right */
    for (j = 1; j < (img_dim.np+1); j++) {
      old[0][j]    = old[img_dim.mp][j];
      old[img_dim.mp+1][j] = old[1][j];
    }

    /* reconstruct image, halo swap not needed in serial */
    for (i = 1; i < (img_dim.mp+1); i++) {
        for (j = 1; j < (img_dim.np+1); j++) {
            new[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
        }
    }

    /* set old = new for next iteration, while finding the max delta value */
    for (i = 1; i < (img_dim.mp + 1); i++) {
        for (j = 1; j < (img_dim.np + 1); j++) {
            delta = fabs(new[i][j] - old[i][j]);
            if (delta > retval.delta) {
                retval.delta = delta;
            }
            retval.sum += new[i][j];
            old[i][j] = new[i][j];
        }
    }
    return retval;
}
