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
 * @file serial/comms.c
 * @author James Clark
 * @brief Serial Communication Code
 */

#include <math.h>
#include <mpi.h>

#include <precision.h>
#include <functions.h>

void get_cart_comm (int * rank, int * size, int * dims, MPI_Comm * cart_comm) {
    *cart_comm = (MPI_Comm) 0;
    dims[0] = 1;
    dims[1] = 1;
}

/* set local to global for serial */
void scatter_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global) {
    int i, j;
    for (i = 0; i < img_dim.mp; i++) {
        for (j = 0; j < img_dim.np; ++j) {
            local[i+1][j+1] = global[i][j];
        }
    }
}

void gather_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global) {
    int i, j;
    for (i = 0; i < img_dim.mp; i++) {
        for (j = 0; j < img_dim.np; ++j) {
            global[i][j] = local[i+1][j+1];
        }
    }
}

/* No reduce is needed in serial, just give back what was given */
void reduce (MPI_Comm cart_comm, MPI_Op op, real * local, real * global) {
    *global = *local;
}
