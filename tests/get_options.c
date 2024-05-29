/**
 *
 * @file get_options.c
 *
 * @copyright 2006-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @date 2023-12-11
 *
 */
#include <spm_tests.h>
#include <unistd.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

/**
 * @brief Print default usage for SpM binaries
 */
static inline void
spm_usage(void)
{
    fprintf(stderr,
            "Matrix input (mandatory):\n"
            " -0 --rsa     : RSA/Matrix Market Fortran driver (only real)\n"
            " -1 --hb      : Harwell Boeing C driver\n"
            " -2 --ijv     : IJV coordinate C driver\n"
            " -3 --mm      : Matrix Market C driver\n"
            " -4 --spm     : SPM Matrix driver\n"
            " -9 --lap     : Generate a Laplacian (5-points stencil)\n"
            " -x --xlap    : Generate an extended Laplacian (9-points stencil)\n"
            " -G --graph   : SCOTCH Graph file\n"
            "\n"
            "Matrix input (optional):\n"
            " -a --scatter : Use the distributed generator or scatter the input\n"
            "                spm when SpM is compiled with MPI (default: 0)\n"
            "                0: Replicate the spm, 1: Scatter the spm\n"
            " -d --doftype : cX -> Set to constant dof of size X\n"
            "                vX -> Set to random variadic dof in the range [1, X]\n"
            "                If -d is not specified, single dof are used (equiv. to c1 or v1)\n"
            "\n"
            " -h --help    : this message\n"
            "\n"
            );
}

/**
 * @brief Define the options and their requirement used by SpM
 */
#define GETOPT_STRING "0:1:2:3:4:9:x:G:a:d:h"

#if defined(HAVE_GETOPT_LONG)
/**
 * @brief Define the long options when getopt_long is available
 */
static struct option long_options[] =
{
    {"rsa",         required_argument,  0, '0'},
    {"hb",          required_argument,  0, '1'},
    {"ijv",         required_argument,  0, '2'},
    {"mm",          required_argument,  0, '3'},
    {"spm",         required_argument,  0, '4'},
    {"lap",         required_argument,  0, '9'},
    {"xlap",        required_argument,  0, 'x'},
    {"graph",       required_argument,  0, 'G'},

    {"scatter",     required_argument,  0, 'a'},
    {"doftype",     required_argument,  0, 'd'},
    {"help",        no_argument,        0, 'h'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

/**
 *******************************************************************************
 *
 * @ingroup spm_examples
 *
 * @brief SpM helper function to read command line options in examples.
 *
 * This function takes the command line arguments, and read the given parameters
 * (integers and doubles), as well as the matrix filename and the driver to read
 * it.
 *
 *******************************************************************************
 *
 * @param[in] argc
 *          The number of input parameters
 *
 * @param[in] argv
 *          The NULL terminated list of parameters
 *
 * @param[out] options
 *          The test options filled by the reading of the parameters (driver,
 *          filename, dofs... )
 *
 *******************************************************************************/
void
spmGetOptions( int          argc,
               char       **argv,
               spm_test_t  *options )
{
    int c;

    if (argc == 1) {
        spm_usage();
        exit(0);
    }

    options->driver  = (spm_driver_t)-1;
    options->doftype = 'c';
    options->dofmax  = 1;
    options->spmdist = 0;
    do
    {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long( argc, argv, GETOPT_STRING,
                         long_options, NULL );
#else
        c = getopt( argc, argv, GETOPT_STRING );
#endif  /* defined(HAVE_GETOPT_LONG) */

        switch(c)
        {
        case '0':
#if defined(SPM_WITH_FORTRAN)
            options->driver   = SpmDriverRSA;
            options->filename = strdup( optarg );
#else
            fprintf(stderr, "spmGetOptions: Please compile with SPM_WITH_FORTRAN option to enable RSA driver or use HB driver instead\n");
            goto unknown_option;
#endif
            break;

        case '1':
            options->driver   = SpmDriverHB;
            options->filename = strdup( optarg );
            break;

        case '2':
            options->driver   = SpmDriverIJV;
            options->filename = strdup( optarg );
            break;

        case '3':
            options->driver   = SpmDriverMM;
            options->filename = strdup( optarg );
            break;

        case '4':
            options->driver   = SpmDriverSPM;
            options->filename = strdup( optarg );
            break;

        case '9':
            options->driver   = SpmDriverLaplacian;
            options->filename = strdup( optarg );
            break;

        case 'x':
            options->driver   = SpmDriverXLaplacian;
            options->filename = strdup( optarg );
            break;

        case 'G':
            options->driver   = SpmDriverGraph;
            options->filename = strdup( optarg );
            break;

        case 'd':
            if ( sscanf( optarg, "%c%d", &(options->doftype), &(options->dofmax) ) != 2 ) {
                fprintf(stderr, "\n%s is not a correct value for the dof parameter\n\n", optarg );
                goto unknown_option;
            }
            break;

        case 'a':
            options->spmdist = !( atoi( optarg ) == 0 );
            break;

        case 'h':
            spm_usage();
            exit(EXIT_FAILURE);

        case ':':
            fprintf(stderr, "\nOption %c is missing an argument\n\n", c );
            goto unknown_option;

        case '?': /* getopt_long already printed an error message. */
            spm_usage();
            exit(EXIT_FAILURE);
        default:
            break;
        }
    } while( -1 != c );

    return;

  unknown_option:
    spm_usage();
    exit(EXIT_FAILURE);
}
