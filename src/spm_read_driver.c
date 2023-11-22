/**
 * @file spm_read_driver.c
 *
 * SParse Matrix package file driver.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2023-01-11
 *
 **/
#include "common.h"
#include "spm_drivers.h"
#if defined(SPM_WITH_SCOTCH)
#include <scotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 * @brief Import a sparse matrix from a Scotch data file.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The name of the file that stores the matrix (see driver).
 *
 * @param[inout] spm
 *          On entry, an allocated sparse matrix structure.
 *          On exit, the filled sparse matrix structure with the matrix from the
 *          file.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the file reading happened successfully,
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
static inline int
spm_read_scotch( const char *filename,
                 spmatrix_t *spm )
{
#if !defined(SPM_WITH_SCOTCH)
    (void)filename;
    (void)spm;

    fprintf( stderr, "Scotch driver to read graph file unavailable.\n"
             "Compile with Scotch support to provide it\n" );
    return SPM_ERR_BADPARAMETER;

#else

    SCOTCH_Graph sgraph;
    FILE        *file;
    SCOTCH_Num   baseval = 1;

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

    SCOTCH_graphLoad( &sgraph, file, -1, 0 );
    SCOTCH_graphData( &sgraph, &baseval, &(spm->n), &(spm->colptr), NULL, NULL, NULL,
                      &(spm->nnz), &(spm->rowptr), NULL );
    fclose(file);

    spm->baseval = baseval;
    spm->mtxtype = SpmGeneral;
    spm->flttype = SpmPattern;
    spm->fmttype = SpmCSC;
    spm->dof     = 1;

    spmUpdateComputedFields( spm );
    return SPM_SUCCESS;
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Import a matrix file into a spm structure for a specific communicator.
 *
 * This function read or generate a sparse matrix from a file to store it into a
 * spm structure. The different formats accepted by this driver are described by
 * the driver field.
 *
 *******************************************************************************
 *
 * @param[in] scatter
 *          Boolean to specify if the final spm must be scattered or not.
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
 *          The MPI communicator of the problem
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the file reading happened successfully,
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
static inline int
spm_read_driver( int          scatter,
                 spm_driver_t driver,
                 const char  *filename,
                 spmatrix_t  *spm,
                 SPM_Comm     comm )
{
    int is_centralized = 1;
    int rc = SPM_SUCCESS;

    if ( filename == NULL ) {
        fprintf( stderr, "spmReadDriver[Dist]: invalid filename parameter\n" );
        return SPM_ERR_BADPARAMETER;
    }

    if ( spm == NULL ) {
        fprintf( stderr, "spmReadDriver[Dist]: invalide spm parameter\n" );
        return SPM_ERR_BADPARAMETER;
    }

    spmInitDist( spm, comm );

    switch(driver)
    {
    case SpmDriverRSA:
        /* The RSA driver is no longer supported in fortran */
        fprintf(stderr, "RSA driver is no longer supported and is replaced by the HB driver\n");
        spm_attr_fallthrough;

    case SpmDriverHB:
        /* TODO: Possible to read the RHS, the solution or a guess of the solution */
        rc = readHB( filename, spm );
        break;

    case SpmDriverIJV:
        rc = readIJV( filename, spm );
        break;

    case SpmDriverMM:
        rc = readMM( filename, spm );
        break;

    case SpmDriverLaplacian:
        rc = genLaplacian( filename, spm );
        is_centralized = 0;
        break;

    case SpmDriverXLaplacian:
        rc = genExtendedLaplacian( filename, spm );
        is_centralized = 0;
        break;

    case SpmDriverSPM:
        rc = spmLoad( spm, filename );
        break;

    case SpmDriverGraph:
        rc = spm_read_scotch( filename, spm );
        break;

    default:
        fprintf(stderr, "spmReadDriver: Driver not implemented\n");
        return SPM_ERR_UNKNOWN;
    }

#if defined(SPM_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                   MPI_MAX, comm );
#endif
    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr,"spmReadDriver[Dist]: error while reading the input %s\n", filename );
        return rc;
    }

#if defined(SPM_WITH_MPI)
    if ( spm->clustnbr > 1 ) {

        if ( is_centralized && scatter )
        {
            /* Scatter the spm among the processes */
            spmatrix_t spm_dist;

            spmScatter( &spm_dist, -1, spm, 0, NULL, 1, spm->comm );

            /* Switch the data structure */
            spmExit( spm );
            memcpy( spm, &spm_dist, sizeof(spmatrix_t) );
        }

        if ( !is_centralized && !scatter )
        {
            /* Gather the spm to replicate it on each node */
            spmatrix_t spm_glob;

            spmGather( spm, -1, &spm_glob );

            /* Switch the data structure */
            spmExit( spm );
            memcpy( spm, &spm_glob, sizeof(spmatrix_t) );
        }
    }
#endif

    (void)is_centralized;
    (void)scatter;
    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Import a matrix file into an spm structure for a specific communicator.
 *
 * This function read or generate a sparse matrix from a file to store it into
 * an spm structure. The different formats accepted by this driver are described
 * by the driver field.
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
 *          The MPI communicator of the problem
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the file reading happened successfully,
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
int
spmReadDriverDist( spm_driver_t driver,
                   const char  *filename,
                   spmatrix_t  *spm,
                   SPM_Comm     comm )
{
    return spm_read_driver( 1, driver, filename,
                            spm, comm );
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Import a matrix file into a spm structure.
 *
 * This function read or generate a sparse matrix from a file to store it into
 * an spm structure. The different formats accepted by this driver are described
 * by the driver field.
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
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the file reading happened successfully,
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
int
spmReadDriver( spm_driver_t driver,
               const char  *filename,
               spmatrix_t  *spm )
{
    return spm_read_driver( 0, driver, filename,
                            spm, MPI_COMM_WORLD );
}
