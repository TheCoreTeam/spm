/**
 * @file spm_read_driver.c
 *
 * SParse Matrix package file driver.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "spm_drivers.h"
#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup spm_spm
 *
 * @brief Import a matrix file into a spm structure.
 *
 * This function read or generate a sparse matrix from a file to store it into a
 * spm structure. The different formats accepted by this driver are described by
 * the driver field.
 *
 *******************************************************************************
 *
 * @param[in] driver
 *          This defines the driver to use to create the spm structure:
 *          - SpmDriverRSA
 *          - SpmDriverHB
 *          - SpmDriverIJV
 *          - SpmDriverMM
 *          - SpmDriverLaplacian
 *          - SpmDriverXLaplacian
 *          - SpmDriverGraph
 *          - SpmDriverSPM
 *
 * @param[in] filename
 *          The name of the file that stores the matrix (see driver).
 *
 * @param[inout] spm
 *          On entry, an allocated sparse matrix structure.
 *          On exit, the filled sparse matrix structure with the matrix from the
 *          file.
 *
 * @param[in] comm
 *          The MPI communicator on which to distribute the sparse matrix. This
 *          is also used in case of distributed formats.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the file reading happened successfully,
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
int
spmReadDriver( spm_driver_t  driver,
               const char   *filename,
               spmatrix_t   *spm,
               MPI_Comm      comm )
{
    int mpirank = 0;
    /* int mpiinit; */

    spmInit(spm);

    /* MPI_Initialized( &mpiinit ); */
    /* if (mpiinit) { */
    /*     MPI_Comm_rank( comm, &mpirank ); */
    /* } */
    mpirank = 0;

    if ( mpirank == 0 )
    {
        switch(driver)
        {
#if defined(SPM_WITH_FORTRAN)
        case SpmDriverRSA:
            readRSA( filename, spm );
            break;
#endif
        case SpmDriverHB:
            /* TODO: Possible to read the RHS, the solution or a guess of the solution */
            readHB( filename, spm );
            break;

        case SpmDriverIJV:
            readIJV( filename, spm );
            break;

        case SpmDriverMM:
            readMM( filename, spm );
            break;

        case SpmDriverLaplacian:
            genLaplacian( filename, spm );
            break;

        case SpmDriverXLaplacian:
            genExtendedLaplacian( filename, spm );
            break;

        case SpmDriverSPM:
        {
            FILE *file = fopen( filename, "r" );
            if ( file == NULL ) {
                fprintf( stderr,"spmReadDriver: impossible to open the file %s\n", filename );
                return SPM_ERR_FILE;
            }

            if ( spmLoad( spm, file ) != SPM_SUCCESS ) {
                fprintf( stderr,"spmReadDriver: error while reading the file %s\n", filename );
                return SPM_ERR_FILE;
            }

            fclose( file );
        }
        break;

        case SpmDriverGraph:
#if defined(HAVE_SCOTCH)
        {
            SCOTCH_Graph sgraph;
            FILE *file;

            file = fopen( filename, "r" );
            if ( file == NULL ) {
                fprintf( stderr,"spmReadDriver: impossible to open the file %s\n", filename );
                return SPM_ERR_FILE;
            }

            /* Check integer compatibility */
            if (sizeof(spm_int_t) != sizeof(SCOTCH_Num)) {
                fprintf( stderr,"Inconsistent integer type\n");
                fclose(file);
                return SPM_ERR_INTEGER_TYPE;
            }

            SCOTCH_graphLoad( &sgraph, file, 1, 0 );
            SCOTCH_graphData( &sgraph, NULL, &(spm->n), &(spm->colptr), NULL, NULL, NULL,
                              &(spm->nnz), &(spm->rowptr), NULL );
            fclose(file);

            spm->mtxtype = SpmGeneral;
            spm->flttype = SpmPattern;
            spm->fmttype = SpmCSC;
            spm->dof = 1;
            spmUpdateComputedFields( spm );
        }
#else
        {
            fprintf(stderr, "Scotch driver to read graph file unavailable.\n"
                    "Compile with Scotch support to provide it\n");
            return SPM_ERR_BADPARAMETER;
        }
#endif
        break;

        default:
            fprintf(stderr, "spmReadDriver: Driver not implemented\n");
        }
    }

    /* if ( mpiinit ) */
    /* { */
    /*     spm_int_t nnz = spm->nnz; */

    /*     MPI_Bcast( spm, sizeof(spmatrix_t), MPI_CHAR, 0, comm ); */
    /*     MPI_Bcast( &nnz, 1, SPM_MPI_INT, 0, comm ); */

    /*     if ( mpirank != 0 ) */
    /*     { */
    /*         spm->colptr = (spm_int_t *) malloc((spm->gN+1) * sizeof(spm_int_t)); */
    /*         spm->rowptr = (spm_int_t *) malloc(nnz * sizeof(spm_int_t)); */
    /*         spm->loc2glob = NULL; */
    /*         spm->loc2glob = NULL; */
    /*     } */

    /*     MPI_Bcast( spm->colptr, spm->gN+1, SPM_MPI_INT, 0, comm ); */
    /*     MPI_Bcast( spm->rowptr, nnz,       SPM_MPI_INT, 0, comm ); */

    /*     if ( spm->flttype != SpmPattern ) { */
    /*         size_t eltsize = spm_size_of( spm->flttype ); */
    /*         if ( mpirank != 0 ) { */
    /*             spm->values = (void *) malloc(nnz * eltsize); */
    /*         } */
    /*         MPI_Bcast( spm->values, nnz * eltsize, MPI_CHAR, 0, comm ); */
    /*     } */
    /* } */

    spmUpdateComputedFields( spm );

    (void)comm;
    return SPM_SUCCESS;
}
