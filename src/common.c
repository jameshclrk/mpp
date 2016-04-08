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
 * @file common.c
 * @author James Clark
 * @brief Code Common to both parallel and serial
 */

#include <stdio.h>
#include <mpi.h>

#include <pgmio.h>
#include <precision.h>
#include <functions.h>

/**
 * @brief Calculate the boundary values for the sawtooth
 * @param i the current location in the image
 * @param m the width of the image
 * @return return the value for the current point on the boundary
 */
real boundaryval (int i, int m) {
  real val;

  val = 2.0*((real)(i-1))/((real)(m-1));
  if (i >= m/2+1) val = 2.0-val;

  return val;
}

/**
 * @brief Set up the initial guess and sawtooth boundary
 * @param cart_comm the cartesian communicator for the processes
 * @param rank the rank of the process calling the function
 * @param img_dim the dimensions of the local and global data
 * @param old array to make the initial guess in and to calculate sawtooth
 */
void setup_reconstruct (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** old){
    int i, j;

    /* set initial guess to white (255), including halos */
    for (i = 0; i < (img_dim.mp + 2); i++) {
        for (j = 0; j < (img_dim.np + 2); j++) {
            old[i][j] = 255.0;
        }
    }
    /* calculate the sawtooth halos */
    sawtooth(cart_comm, rank, img_dim, old);
}

/**
 * @brief Get the dimensions of an image (Wrapper for pgmsize)
 * @param filename the file to find the dimensions of
 * @param nx pointer to store the x dimension
 * @param ny pointer to store the y dimension
 */
void image_size (char *filename, int *nx, int *ny) {
    pgmsize(filename, nx, ny);
}

/**
 * @brief Read in a PGM file to an array. Only rank 0 can read (Wrapper for pgmread)
 * @param rank the rank of the calling process
 * @param filename the file to read
 * @param img_dim the dimensions of the image
 * @param data the array to read the image in to
 */
void image_read (int rank, char * filename, image_dimensions img_dim, real ** data) {
    if (rank == 0) {
        printf("Reading %d x %d picture from file: %s\n", img_dim.m, img_dim.n, filename);
        pgmread(filename, &data[0][0], img_dim.m, img_dim.n);
    }
}

/**
 * @brief Write a PGM file from an array. Only rank 0 can write (Wrapper for pgmwrite)
 * @param rank the rank of the calling process
 * @param filename the file to write to
 * @param img_dim the dimensions of the image
 * @param data the array to write to disk
 */
 void image_write (int rank, char * filename, image_dimensions img_dim, real ** data) {
    if (rank == 0) {
        pgmwrite(filename, &data[0][0], img_dim.m, img_dim.n);
    }

}
