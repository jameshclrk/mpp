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
 * @file parallel/comms.c
 * @author James Clark
 * @brief Parallel Communication Code
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include <precision.h>
#include <functions.h>

/**
 * @brief Create the cartesian communicator and dimensions for the topology.
 * @param rank the rank of the process calling the function
 * @param size the number of processes in the communicator
 * @param dims stores how many processes are in each dimension
 * @param cart_comm the cartesian communicator for the processes
 */
void get_cart_comm (int * rank, int * size, int * dims, MPI_Comm * cart_comm) {
    /* periodic only in one dimension */
    int periods[2]  = {1,0};
    /* let MPI decide the best dimensions */
    MPI_Dims_create(*size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, cart_comm);
    /* double check rank, in case Cart_create reordered them */
    MPI_Comm_rank(*cart_comm, rank);
}

/**
 * @brief Scatters global, from process 0, to local, on all proceses.
 * @param cart_comm the cartesian communicator for the processes
 * @param rank the rank of the process calling the function
 * @param size the number of processes in the communicator
 * @param img_dim the dimensions of the local and global data
 * @param local where the data is being scattered to
 * @param global the data source
 */
void scatter_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global) {
    /* Send data back to rank 0 */
    int i;
    int offset_m, offset_n;
    int num_req = 0;
    int coords[2] = {0,0};
    MPI_Request requests[size+1];
    MPI_Status  statuses[size+1];
    MPI_Datatype send_array_type;
    MPI_Datatype recv_array_type;

    /* derived type for sending and receiving */
    MPI_Type_vector(img_dim.mp, img_dim.np, img_dim.n,  MPI_REALNUM, &send_array_type);
    MPI_Type_vector(img_dim.mp, img_dim.np, (img_dim.np)+2,  MPI_REALNUM, &recv_array_type);
    MPI_Type_commit(&send_array_type);
    MPI_Type_commit(&recv_array_type);

    if (rank != 0 ) {
        /* all ranks (other than zero) receive from zero */
        MPI_Irecv(&local[1][1], 1, recv_array_type, 0, rank, cart_comm, &requests[0]);
        num_req++;
    } else {
        /* rank zero receives from itself, then sends to everyone */
        MPI_Irecv(&local[1][1], 1, recv_array_type, 0, rank, cart_comm, &requests[size]);
        num_req++;
        for (i = 0; i < size; i++) {
            /* get the coords of the receiver */
            MPI_Cart_coords(cart_comm, i, 2, coords);
            /* calculate the where the data should be sent from */
            offset_m = coords[0]*img_dim.mp;
            offset_n = coords[1]*img_dim.np;
            MPI_Issend(&global[offset_m][offset_n], 1, send_array_type, i, i, cart_comm, &requests[i]);
            num_req++;
        }
    }
    MPI_Waitall(num_req, requests, statuses);
}

/**
 * @brief Gathers local from all proceses to global on process 0.
 * @param cart_comm the cartesian communicator for the processes
 * @param rank the rank of the process calling the function
 * @param size the number of processes in the communicator
 * @param img_dim the dimensions of the local and global data
 * @param local where the data is being scattered to
 * @param global the data source
 */
void gather_data (MPI_Comm cart_comm, int rank, int size, image_dimensions img_dim, real ** local, real ** global) {
    int i;
    int offset_m, offset_n;
    int num_req = 0;
    int coords[2] = {0,0};
    MPI_Request requests[size+1];
    MPI_Status  statuses[size+1];
    MPI_Datatype send_array_type;
    MPI_Datatype recv_array_type;

    /* derived type for sending and receiving */
    MPI_Type_vector(img_dim.mp, img_dim.np, img_dim.np+2,  MPI_REALNUM, &send_array_type);
    MPI_Type_vector(img_dim.mp, img_dim.np, img_dim.n,   MPI_REALNUM, &recv_array_type);
    MPI_Type_commit(&send_array_type);
    MPI_Type_commit(&recv_array_type);

    if (rank != 0 ) {
        /* all ranks (other than zero) sent to zero */
        MPI_Issend(&local[1][1], 1, send_array_type, 0, rank, cart_comm, &requests[0]);
        num_req++;
    } else {
        /* rank zero sends to itself, then receives from everyone */
        MPI_Issend(&local[1][1], 1, send_array_type, 0, rank, cart_comm, &requests[size]);
        num_req++;
        for (i = 0; i < size; i++) {
            /* get the coords of the sender */
            MPI_Cart_coords(cart_comm, i, 2, coords);
            /* calculate the where the data should be received to */
            offset_m = coords[0]*img_dim.mp;
            offset_n = coords[1]*img_dim.np;
            MPI_Irecv(&global[offset_m][offset_n], 1, recv_array_type, i, i, cart_comm, &requests[i]);
            num_req++;
        }
    }
    MPI_Waitall(num_req, requests, statuses);
}

/**
 * @brief Performs an Allreduce for a real number.
 * @param cart_comm the cartesian communicator for the processes
 * @param op the operation to perform on the reduce
 * @param local each processes local number
 * @param global the number to store the reduced version
 */
void reduce (MPI_Comm cart_comm, MPI_Op op, real * local, real * global) {
    MPI_Allreduce(local, global, 1, MPI_REALNUM, op, cart_comm);
}
