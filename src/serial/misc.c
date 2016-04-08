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
 * @file serial/misc.c
 * @author James Clark
 * @brief Misc Serial Code
 */

#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

#include <precision.h>
#include <functions.h>

void sawtooth (MPI_Comm cart_comm, int rank, image_dimensions img_dim, real ** old) {
    int i;
    real val;
    for (i = 1; i < (img_dim.mp + 1); i++) {
      /* compute sawtooth value */
      val = boundaryval(i, img_dim.m);

      old[i][0]   = 255.0*val;
      old[i][img_dim.np+1] = 255.0*(1.0-val);
    }
}

double get_time () {
    struct timeval mtime;
    gettimeofday(&mtime, NULL);
    return ((mtime.tv_sec) * 1000000u + mtime.tv_usec) / 1.e6;
}

/* set up rank/size */
void init (int argc, char * argv[], int * rank, int * size) {
    *rank = 0;
    *size = -1;
}

void finalise () {
    return;
}

void m_abort () {
    exit(1);
}
