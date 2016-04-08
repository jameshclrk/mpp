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


/* The functions and variables included here conform to the argp API */

const char *argp_program_version =
  "reconstruct 0.1";

/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    args *arguments = state->input;

    switch (key)
    {
        case 'i':
            arguments->iterations = atoi(arg);
            break;
        case 's':
            arguments->step = atoi(arg);
            break;
        case 'd':
            arguments->delta = atof(arg);
            break;
        case 'o':
            arguments->output = arg;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 1)
            {
                argp_usage(state);
            }
            arguments->filename = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1)
            {
                argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
/* Description of the required arguments */
static char args_doc[] = "edge_file";

/* Description of optional arguments
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
 */
static struct argp_option options[] =
{
  {"max_iterations", 'i', "MAX_ITERATIONS", 0, "Manually specify the maximum number of iterations"},
  {"step", 's', "STEP", 0, "Manually specify the output step"},
  {"delta", 'd', "DELTA", 0, "Manually specify the minimum delta value"},
  {"output_file", 'o', "FILE", 0, "Manually specify the output file"},
  {0}
};
/* Documentation String */
static char doc[] = "reconstruct an image from an edge file.\vDeveloped by:\tB081060";

/*   The ARGP structure itself. */
static struct argp argp = {options, parse_opt, args_doc, doc};
