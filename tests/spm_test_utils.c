/**
 *
 * @file spm_test_utils.c
 *
 * Utils routines to factorize the test files
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-12-11
 *
 **/
#include "spm_tests.h"

const char *fltnames[]   = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
const char *fmtnames[]   = { "CSC", "CSR", "IJV" };
const char *mtxnames[]   = { "General", "Symmetric", "Hermitian" };
const char *dofnames[]   = { "None", "Constant", "Variadic" };
const char *transnames[] = { "NoTrans", "Trans", "ConjTrans" };

/**
 *******************************************************************************
 *
 * @brief This routine will get the options from the command line and initialize
 *        the spm structure
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, an allocated sparse matrix structure.
 *          On exit, the filled sparse matrix structure with the matrix from the
 *          file.
 *
 * @param[in] argc
 *          The number of input parameters
 *
 * @param[in] argv
 *          The NULL terminated list of parameters
 *
 *******************************************************************************
 *
 * @retval the return value of spmReadDriver.
 *
 *******************************************************************************/
int
spmTestGetSpm( spmatrix_t *spm,
               int         argc,
               char      **argv )
{
    spm_test_t options;
    int        rc;

    spmGetOptions( argc, argv, &options );

    if ( options.driver == (spm_driver_t)-1 ) {
        fprintf( stderr, "[%s] Incorrect driver type. Please specify a correct driver.\n", __func__ );
        return SPM_ERR_BADPARAMETER;
    }

    if ( options.spmdist == 1 ) {
        rc = spmReadDriverDist( options.driver, options.filename, spm, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( options.driver, options.filename, spm );
    }
    free(options.filename);
    if ( rc != SPM_SUCCESS ) {
        return rc;
    }

    if ( options.dofmax > 1 ) {
        spmatrix_t spm2;
        int type   = (options.doftype == 'v');
        int dofmax = options.dofmax;

        rc = spmDofExtend( spm, type, dofmax, &spm2 );

        spmExit( spm );
        memcpy( spm, &spm2, sizeof(spmatrix_t) );
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Indicates if we have to skip the loop iteration in the test.
 *
 *******************************************************************************
 *
 * @param[in] flttype
 *          The floating type of the spm values.
 *
 * @param[in] spm_mtxtype
 *          The previous matrix type of the spm.
 *
 * @param[in] new_mtxtype
 *          The next matrix type of the spm in the for loop.
 *
 *******************************************************************************
 *
 * @retval 1 if we need to continue in the for loop, 0 otherwise.
 *
 *******************************************************************************/
int
spmTestPassMtxtype( spm_coeftype_t flttype,
                    spm_mtxtype_t  spm_mtxtype,
                    spm_mtxtype_t  new_mtxtype )
{
    int pass       = 0;
    int is_complex = ( (flttype == SpmComplex64) || (flttype == SpmComplex32) );

    if ( (new_mtxtype == SpmHermitian) &&
        ((spm_mtxtype != SpmHermitian) || !is_complex) )
    {
        pass = 1;
    }

    if ( (new_mtxtype != SpmGeneral) &&
         (spm_mtxtype == SpmGeneral) )
    {
        pass = 1;
    }

    return pass;
}

/**
 *******************************************************************************
 *
 * @brief End of each test.
 *
 *******************************************************************************
 *
 * @param[in] err
 *          The amount of errors in the test.
 *
 * @param[in] clustnum
 *          cluster number of the spm.
 *
 *******************************************************************************
 *
 * @retval EXIT_SUCCESS if err = 0, EXIT_FAILURE otherwise.
 *
 *******************************************************************************/
int
spmTestEnd( int err,
            int clustnum )
{
    if ( err == 0 ) {
        if ( clustnum == 0 ) {
            printf(" -- All tests PASSED --\n");
        }
        return EXIT_SUCCESS;
    }
    else {
        printf("[%d] -- %d tests FAILED --\n", clustnum, err);
        return EXIT_FAILURE;
    }
}

/**
 *******************************************************************************
 *
 * @brief Create a loc2glob array
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm that will allow us to create the loc2glob.
 *
 * @param[out] loc2globptr
 *          The pointer to the loc2glob array.
 *
 * @param[in] l2gtype
 *          Type of the loc2glob ( SpmContinuous, SpmRoundRobin )
 *
 *******************************************************************************
 *
 * @retval The size of the loc2glob array.
 *
 *******************************************************************************/
spm_int_t
spmTestCreateL2g( const spmatrix_t *spm,
                  spm_int_t       **loc2globptr,
                  spm_l2gtype_t     l2gtype )
{
    spm_int_t *loc2glob;
    spm_int_t  i, baseval, size;
    int        clustnum, clustnbr;

    baseval  = spm->baseval;
    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;
    size     = spm->gN / clustnbr;
    if ( l2gtype == SpmRoundRobin ) {
        spm_int_t ig = clustnum;
        if ( clustnum < (spm->gN % clustnbr) ) {
            size++;
        }
        loc2glob = malloc( size * sizeof(spm_int_t) );
        *loc2globptr = loc2glob;

        for ( i=0; i<size; i++, loc2glob++, ig+=clustnbr )
        {
            *loc2glob = ig + baseval;
        }
    }
    else { /* SpmContinuous or SpmRandom */
        size = spm_create_loc2glob_continuous( spm, loc2globptr );
    }

    if ( l2gtype == SpmRandom ) {
        spm_int_t idx1, idx2, tmp;
        i        = 0;
        loc2glob = *loc2globptr;
        while ( i < ( size / 2 ) )
        {
            /* Get 2 randoms index */
            idx1 = rand() % size;
            idx2 = rand() % size;

            /* Swap values */
            tmp            = loc2glob[idx1];
            loc2glob[idx1] = loc2glob[idx2];
            loc2glob[idx2] = tmp;

            i++;
        }
    }

    return size;
}

/**
 *******************************************************************************
 *
 * @brief Convert the spm and dump the matrix in a generated file
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm that will be converted.
 *
 * @param[in] newtype
 *          The new type for the spm
 *
 * @param[in] cycle
 *          String that indicates the current cycle (cycle1, cycle2 or end)
 *
 *******************************************************************************
 *
 * @retval The error code of the spmConvert routine
 *
 *******************************************************************************/
int
spmTestConvertAndPrint( spmatrix_t   *spm,
                        spm_fmttype_t newtype,
                        const char   *cycle )
{
    char *filename;
    FILE *f;
    int   rc, err, baseval;

    printf( "   -- Test Conversion %s -> %s: ", fmtnames[ spm->fmttype - SpmCSC ], fmtnames[ newtype - SpmCSC ] );
    rc  = spmConvert( newtype, spm );
    err = ( rc != SPM_SUCCESS ) || ( spm->fmttype != newtype );
    baseval = spm->baseval;
    if ( spm->loc2glob != NULL ) {
        rc = asprintf( &filename,
                       "convert_dist_b%d_%s_%s_%s_%d.dat",
                       baseval,
                       mtxnames[spm->mtxtype - SpmGeneral],
                       fmtnames[newtype - SpmCSC],
                       cycle,
                       spm->clustnum );
    }
    else {
        rc = asprintf( &filename, "convert_b%d_%s_%s_%s.dat",
                       baseval,
                       mtxnames[spm->mtxtype - SpmGeneral],
                       fmtnames[newtype - SpmCSC],
                       cycle );
    }

    if ( ( f = fopen( filename, "w" ) ) == NULL ) {
        fprintf( stderr, "spm_convert_test:%s:%s", cycle, fmtnames[newtype - SpmCSC] );
        free( filename );
        return EXIT_FAILURE;
    }
    spmPrint( spm, f );
    fclose( f );
    free( filename );

    return err;
}

/**
 *******************************************************************************
 *
 * @brief Realize the main loop of the test that concerns one matrix.
 *
 *******************************************************************************
 *
 * @param[inout] original
 *          The original spm to test.
 *
 * @param[in] spm_test_check
 *          Routine that will test the spm.
 *
 * @param[in] to_scatter
 *          Boolean that indicates if we have to scatter the spm.
 *
 *******************************************************************************
 *
 * @retval The amount of errors in the loop.
 *
 *******************************************************************************/
int
spmTestLoop( spmatrix_t        *original,
             spm_test_check_fct spm_test_check,
             int                to_scatter )
{
    spmatrix_t   *spm, spmtmp;
    spm_fmttype_t fmttype;
    spm_mtxtype_t origtype, mtxtype;
    int           rc, err=0;
    int           clustnum, baseval;
    int           dofidx;

    if ( original->dof == 1 ) {
        dofidx = 0;
    }
    else if ( original->dof > 1 ) {
        dofidx = 1;
    }
    else {
        dofidx = 2;
    }

    origtype = original->mtxtype;
    clustnum = original->clustnum;
    /* Loop on the fmttype : SpmCSC, SpmCSR, SpmIJV*/
    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
    {
        if ( spmConvert( fmttype, original ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to convert to %s format\n", fmtnames[fmttype] );
            continue;
        }

        /* Scatter the Spm if we are on a distributed test */
        if ( to_scatter ) {
            int distbycol = (fmttype == SpmCSR) ? 0 : 1;

            rc = spmScatter( &spmtmp, -1, original, -1, NULL, distbycol, MPI_COMM_WORLD );
            if ( rc != SPM_SUCCESS ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }
            spm = &spmtmp;
        }
        else {
            spm = original;
        }

        /*
         * If the format is IJV, let's generate once and for all the glob2loc
         * used by z_spmRhsGenRndDist
         */
        if ( spm->fmttype == SpmIJV ) {
            spm_getandset_glob2loc( spm );
        }

        /* Loop on the baseval */
        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( spm, baseval );

            for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
            {
                if ( spmTestPassMtxtype( spm->flttype, origtype, mtxtype ) ) {
                    continue;
                }
                spm->mtxtype = mtxtype;

                if( clustnum == 0 ) {
                    printf(" Case: %s / %s / %d / %s / %s ",
                            fltnames[spm->flttype],
                            fmtnames[spm->fmttype],
                            baseval,
                            mtxnames[mtxtype - SpmGeneral],
                            dofnames[dofidx] );
                }

                /* Do the test */
                rc = spm_test_check( spm );
                PRINT_RES( rc );
            }
        }

        if ( to_scatter ) {
            spmExit( spm );
        }
    }

    return err;
}

/**
 *******************************************************************************
 *
 * @brief Realize the main loop of the test that 2 matrices : one centralized,
 *        one distributed.
 *
 *******************************************************************************
 *
 * @param[inout] original
 *          The original spm to scatter and test.
 *
 * @param[in] spm_test_check
 *          Routine that will test and compare the spm.
 *
 *******************************************************************************
 *
 * @retval The amount of errors in the loop.
 *
 *******************************************************************************/
int
spmTestLoop2( spmatrix_t         *original,
              spm_test_check2_fct spm_test_check2 )
{
    spmatrix_t   *spm, spmtmp;
    spm_fmttype_t fmttype;
    spm_mtxtype_t origtype, mtxtype;
    int           rc, err=0;
    int           clustnum, baseval, distbycol;
    int           dofidx;

    if ( original->dof == 1 ) {
        dofidx = 0;
    }
    else if ( original->dof > 1 ) {
        dofidx = 1;
    }
    else {
        dofidx = 2;
    }

    origtype = original->mtxtype;
    clustnum = original->clustnum;
    /* Loop on the fmttype : SpmCSC, SpmCSR, SpmIJV*/
    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
    {
        if ( spmConvert( fmttype, original ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to convert to %s format\n", fmtnames[fmttype] );
            continue;
        }

        /* Scatter the Spm if we are on a distributed test */
        distbycol = (fmttype == SpmCSR) ? 0 : 1;

        if ( original->loc2glob == NULL ) {
            rc = spmScatter( &spmtmp, -1, original, -1, NULL, distbycol, MPI_COMM_WORLD );
            if ( rc != SPM_SUCCESS ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }
        }
        else {
            spmCopy( original, &spmtmp );
        }
        spm = &spmtmp;

        /* Loop on the baseval */
        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( original, baseval );
            spmBase( spm, baseval );

            for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
            {
                if ( spmTestPassMtxtype( spm->flttype, origtype, mtxtype ) ) {
                    continue;
                }
                original->mtxtype = mtxtype;
                spm->mtxtype      = mtxtype;

                if( clustnum == 0 ) {
                    printf(" Case: %s / %s / %d / %s / %s ",
                            fltnames[spm->flttype],
                            fmtnames[spm->fmttype],
                            baseval,
                            mtxnames[mtxtype - SpmGeneral],
                            dofnames[dofidx] );
                }

                /* Do the test */
                rc = spm_test_check2( original, spm );
                PRINT_RES( rc );
            }
        }

        spmExit( spm );
    }

    return err;
}
