/* * MPP Coursework - MPI Edge Reconstruction
 * Copyright (C) 2015,2016 James Clark
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
 * @file parallel/update.c
 * @author James Clark
 * @brief Parallel Update Code
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include <precision.h>
#include <functions.h>

/**
 * @brief Performs one reconstruct operation.
 * @param cart_comm the cartesian communicator for the processes
 * @param rank the rank of the process calling the function
 * @param img_dim the dimensions of the local and global data
 * @param edge stores the original edge data
 * @param old stores the previous operation's data
 * @param new stores the current operation's data
 * @return both the maximum pixel change and the the average pixel value. See ::step_return
 */
step_return update_tick (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** edge, real ** old, real ** new){
    int i, j, i_up, i_down, j_up, j_down;
    real delta = 0.0;
    step_return retval = {0.0, 0.0};

    MPI_Request requests[8];
    MPI_Status  statuses[8];
    MPI_Datatype i_halo;

    /* derived type for halo swaps between horizontal neighbours */
    MPI_Type_vector(img_dim.mp, 1, (img_dim.np)+2,  MPI_REALNUM, &i_halo);
    MPI_Type_commit(&i_halo);

    /* find neighbours */
    MPI_Cart_shift(cart_comm, 0, 1, &j_down, &j_up);
    MPI_Cart_shift(cart_comm, 1, 1, &i_down, &i_up);


    /* non blocking send/recv of halos */
    /* synchronous sends, so data cannot be modifed until send/recv completes */
    MPI_Issend(&old[img_dim.mp][1], img_dim.np, MPI_REALNUM,   j_up, 1, cart_comm, &requests[0]);
    MPI_Issend(&old[1][1],          img_dim.np, MPI_REALNUM, j_down, 2, cart_comm, &requests[1]);
    MPI_Issend(&old[1][img_dim.np],          1,      i_halo,   i_up, 4, cart_comm, &requests[2]);
    MPI_Issend(&old[1][1],                   1,      i_halo, i_down, 3, cart_comm, &requests[3]);

    MPI_Irecv(&old[0][1],            img_dim.np, MPI_REALNUM, j_down, 1, cart_comm, &requests[4]);
    MPI_Irecv(&old[img_dim.mp+1][1], img_dim.np, MPI_REALNUM,   j_up, 2, cart_comm, &requests[5]);
    MPI_Irecv(&old[1][img_dim.np+1],          1,      i_halo,   i_up, 3, cart_comm, &requests[6]);
    MPI_Irecv(&old[1][0],                     1,      i_halo, i_down, 4, cart_comm, &requests[7]);


    /* Rather than waiting for halos, keep doing work by
     * reconstructing the image excluding pixels that need the halos */
    for (i = 2; i < (img_dim.mp); i++) {
        for (j = 2; j < (img_dim.np); j++) {
            new[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
        }
    }

    /* wait for halo swap, hopefully completed by now */
    MPI_Waitall(8, requests, statuses);

    /* reconstruct pixels that depend on halos */
    for (j = 1; j < (img_dim.np + 1); j++) {
        new[1][j] = 0.25 * (old[1][j-1] + old[1][j+1] + old[0][j] + old[2][j] - edge[1][j]);
        new[img_dim.mp][j] = 0.25 * (old[img_dim.mp][j-1] + old[img_dim.mp][j+1] + old[img_dim.mp-1][j] + old[img_dim.mp+1][j] - edge[img_dim.mp][j]);
    }

    for (i = 1; i < (img_dim.mp + 1); i++) {
        new[i][1] = 0.25 * (old[i][0] + old[i][2] + old[i-1][1] + old[i+1][1] - edge[i][1]);
        new[i][img_dim.np] = 0.25 * (old[i][img_dim.np-1] + old[i][img_dim.np+1] + old[i-1][img_dim.np] + old[i+1][img_dim.np] - edge[i][img_dim.np]);
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
