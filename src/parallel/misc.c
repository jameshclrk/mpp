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
 * @file parallel/misc.c
 * @author James Clark
 * @brief Misc Parallel Code
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include <precision.h>
#include <functions.h>

/**
 * @brief Sawtooth boundary condition. Takes the local process's position in the image in to account
 * @param cart_comm the cartesian communicator for the processes
 * @param rank the rank of the process calling the function
 * @param img_dim the dimensions of the local and global data
 * @param old the array to apply the boundary condition to
 */
void sawtooth (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** old) {
    int i;
    int offset_m;
    int coords[2] = {0,0};
    real val;

    /* make sure the coordinates are correct */
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    /* calculate the offset of the local data in the global image */
    offset_m = coords[0]*img_dim.mp;

    /* create the sawtooth value for a local process, taking the global size in to account */
    for (i = 1; i < (img_dim.mp + 1); i++) {
      /* compute sawtooth value */
      val = boundaryval(offset_m + i, img_dim.m);

      old[i][0]   = 255.0*val;
      old[i][img_dim.np+1] = 255.0*(1.0-val);
    }
}

/**
 * @brief Get the current wall time (Wrapper for MPI_Wtime)
 * @return time in seconds since an arbitrary point in the past
 */
double get_time () {
    return MPI_Wtime();
}

/**
 * @brief Initialise MPI
 * @param argc to pass arguments to MPI
 * @param argv to pass arguments to MPI
 * @param rank pointer so the rank can be designated
 * @param size pointer so the size can be designated
 */
void init (int argc, char * argv[], int * rank, int * size) {
    int initialized;

    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, size);
}

/**
 * @brief Finalise MPI
 */
void finalise () {
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
}

/**
 * @brief Abort all MPI processes
 */
void m_abort () {
    MPI_Abort(MPI_COMM_WORLD, 1);
}
