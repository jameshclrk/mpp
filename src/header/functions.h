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
 * @file functions.h
 * @author James Clark
 * @brief Defines the functions and structs.
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

/** Holds the delta and sum of pixels for the update function to return */
typedef struct {
    real delta;  /**< The maximum delta found */
    real sum;    /**< The sum of all pixels */
} step_return;

/** Holds the dimensions for the local and global image data */
typedef struct {
    int m;   /**< The global image size in dim 0 */
    int n;   /**< The global image size in dim 1 */
    int mp;  /**< The local  image size in dim 0 */
    int np;  /**< The local  image size in dim 1 */
} image_dimensions;

/** Holds the arguments for the program */
typedef struct {
    char * filename;      /**< Input file name, required */
    int iterations;       /**< Maximum iterations, provided by -i */
    int step;             /**< Prints information every "step" steps, provided by -s */
    double delta;         /**< Minimum delta value, provided by -d */
    char * output;        /**< Output file name, provided by -o */
} args;

void init (int argc, char * argv[], int * rank, int * size);
void finalise ();
void m_abort ();

double get_time();

step_return update_tick (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** edge, real ** old, real ** new);

void image_size (char *filename, int *nx, int *ny);
void image_read (int rank, char * filename, image_dimensions img_dim, real ** data);
void image_write (int rank, char * filename, image_dimensions img_dim, real ** data);

void get_cart_comm (int * rank, int * size, int * dims, MPI_Comm * cart_comm);
void scatter_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global);
void gather_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global);
void reduce (MPI_Comm cart_comm, MPI_Op op, real * delta, real * global_delta);

void setup_reconstruct (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** old);
void sawtooth (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** old);
real boundaryval (int i, int m);

#endif
