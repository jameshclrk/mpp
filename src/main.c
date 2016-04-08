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
 * @file main.c
 * @author James Clark
 * @brief Contains the main function for the image reconstruct.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <argp.h>
#include <float.h>

#include <arralloc.h>
#include <precision.h>
#include <functions.h>

#include "argp/argp.c"

/* Default options */
/** Defualt output filename */
#define OUTPUT "output.pgm"
/** Default maximum number of iterations */
#define MAX_COUNT 5000
/** Default output interval */
#define STEP 100
/** Default minimum delta value */
#define MIN_DELTA 0.1

int main (int argc, char * argv[]) {
    int rank, size;
    int iteration;
    /* Cartesian dimensions */
    int dims[2] = {0,0};
    /* Struct for global and local image dimensions */
    image_dimensions img_dim;
    /* For timing main loop */
    double t0, t1;
    /* Set inital global values. */
    real global_delta = FLT_MAX,  // Using the max float means the first loop will always occur
         global_average = 1.0;
    /* Pointers for global and local storage */
    real ** main_buf,
         ** edge,
         ** old,
         ** new;
    /* Initialise the return value for the update_step */
    step_return return_val = {1.0, 1.0};

    /* comminucator for cartesian space */
    MPI_Comm cart_comm;

    args arguments;

    /* Set default arguments */
    arguments.iterations = MAX_COUNT;
    arguments.step = STEP;
    arguments.delta = MIN_DELTA;
    arguments.output = OUTPUT;

    /* parse the command line options */
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    /* initialise mpi (if parallel) and get rank and size */
    init(argc, argv, &rank, &size);

    /* Get the cartesian communicator */
    get_cart_comm(&rank, &size, dims, &cart_comm);

    /* confirm to stdout the number of processes and topology */
    if(rank == 0) {
        printf("Running on %d processes\n", size);
        printf("Cartesian topology: %d x %d\n", dims[0], dims[1]);
    }

    /* get the image dimensions */
    image_size(arguments.filename, &(img_dim.m), &(img_dim.n));

    /* check the image can be fit evenly on the processors */
    if ((img_dim.m % dims[0] != 0) || (img_dim.n % dims[1] != 0)) {
        if (rank == 0)
            printf("Cannot fit %dx%d processes evenly on a %dx%d image\n",
                dims[0], dims[1], img_dim.m, img_dim.n);
        m_abort();
    }

    /* get local region dimensions */
    img_dim.mp = img_dim.m/dims[0];
    img_dim.np = img_dim.n/dims[1];

    /* Allocate memory */
    if(rank == 0) {
        /* Only rank 0 needs to allocate the main buffer */
        printf("Allocating memory\n");
        main_buf  = (real **) arralloc(sizeof(real), 2, img_dim.m,    img_dim.n);
    }
    edge      = (real **) arralloc(sizeof(real), 2, img_dim.mp+2, img_dim.np+2);
    old       = (real **) arralloc(sizeof(real), 2, img_dim.mp+2, img_dim.np+2);
    new       = (real **) arralloc(sizeof(real), 2, img_dim.mp+2, img_dim.np+2);

    image_read(rank, arguments.filename, img_dim, main_buf);

    scatter_data(cart_comm, rank, size, img_dim, edge, main_buf);

    setup_reconstruct(cart_comm, rank, img_dim, old);

    if (rank == 0) t0 = get_time();

    /* Reconstruct the image */
    iteration = 0;
    while ((iteration < arguments.iterations) && (global_delta > arguments.delta)) {
        return_val = update_tick(cart_comm, rank, img_dim, edge, old, new);
        reduce(cart_comm, MPI_MAX, &(return_val.delta), &global_delta);

        if (iteration % arguments.step == 0) {
            reduce(cart_comm, MPI_SUM, &(return_val.sum), &global_average);
            if (rank == 0) {
                global_average /=  (img_dim.m * img_dim.n);
                printf("Iteration %7d\tAverage Pixel = %.16f\tGlobal Delta = %.16f\n", iteration, global_average, global_delta);
            }
        }

        iteration++;
    }
    if (rank == 0) {
        t1 = get_time();
        printf("Time for %d iterations: %lf\n", iteration, t1-t0);
    }
    /* reduce global average, just in case the last step was not a stdout step */
    reduce(cart_comm, MPI_SUM, &(return_val.sum), &global_average);
    if (rank == 0) {
        global_average /=  (img_dim.m * img_dim.n);
        printf("Iteration %7d\tAverage Pixel = %.16f\tGlobal Delta = %.16f\n", iteration, global_average, global_delta);
    }

    gather_data(cart_comm, rank, size, img_dim, old, main_buf);

    image_write(rank, arguments.output, img_dim, main_buf);


    /* clean up memory */
    if (rank == 0) free(main_buf);
    free(edge);
    free(old);
    free(new);

    finalise();

    return 0;
}
