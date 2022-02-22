#=

 @file spm.jl

 SPM julia wrapper

 @copyright 2020-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.1.0
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2022-02-22

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_julia

=#

module spm
using CBinding
using Libdl
include("spm_enums.jl")

function spm_library_path()
    x = Libdl.dlext
    return "libspm.$x"
end
libspm = spm_library_path()

@cbindings libspm begin
    @cextern spmInit( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmInitDist( spm::Ptr{spmatrix_t}, comm::__get_mpi_type__() )::Cvoid
end

@cbindings libspm begin
    @cextern spmAlloc( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmExit( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmCopy( spm_in::Ptr{spmatrix_t}, spm_out::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmBase( spm::Ptr{spmatrix_t}, baseval::Cint )::Cvoid
end

@cbindings libspm begin
    @cextern spmFindBase( spm::Ptr{spmatrix_t} )::spm_int_t
end

@cbindings libspm begin
    @cextern spmConvert( ofmttype::Cint, ospm::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmUpdateComputedFields( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmGenFakeValues( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmScatter( spm_scattered::Ptr{spmatrix_t}, root::Cint, opt_spm_gathered::Ptr{spmatrix_t}, opt_n::spm_int_t, opt_loc2glob::Ptr{spm_int_t}, opt_distByColumn::Cint, opt_comm::__get_mpi_type__() )::Cint
end

@cbindings libspm begin
    @cextern spmGather( spm_scattered::Ptr{spmatrix_t}, root::Cint, opt_spm_gathered::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmRedistribute( spm::Ptr{spmatrix_t}, new_n::spm_int_t, newl2g::Ptr{spm_int_t}, newspm::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmNorm( ntype::spm_normtype_t, spm::Ptr{spmatrix_t} )::Cdouble
end

@cbindings libspm begin
    @cextern spmNormVec( ntype::spm_normtype_t, spm::Ptr{spmatrix_t}, x::Ptr{Cvoid}, incx::spm_int_t )::Cdouble
end

@cbindings libspm begin
    @cextern spmNormMat( ntype::spm_normtype_t, spm::Ptr{spmatrix_t}, n::spm_int_t, A::Ptr{Cvoid}, lda::spm_int_t )::Cdouble
end

@cbindings libspm begin
    @cextern spmMatVec( trans::spm_trans_t, alpha::Cdouble, spm::Ptr{spmatrix_t}, x::Ptr{Cvoid}, beta::Cdouble, y::Ptr{Cvoid} )::Cint
end

@cbindings libspm begin
    @cextern spmMatMat( trans::spm_trans_t, n::spm_int_t, alpha::Cdouble, A::Ptr{spmatrix_t}, B::Ptr{Cvoid}, ldb::spm_int_t, beta::Cdouble, C::Ptr{Cvoid}, ldc::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmScal( alpha::Cdouble, spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmScalVec( alpha::Cdouble, spm::Ptr{spmatrix_t}, x::Ptr{Cvoid}, incx::spm_int_t )::Cvoid
end

@cbindings libspm begin
    @cextern spmScalMat( alpha::Cdouble, spm::Ptr{spmatrix_t}, n::spm_int_t, A::Ptr{Cvoid}, lda::spm_int_t )::Cvoid
end

@cbindings libspm begin
    @cextern spmSort( spm::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmMergeDuplicate( spm::Ptr{spmatrix_t} )::spm_int_t
end

@cbindings libspm begin
    @cextern spmSymmetrize( spm::Ptr{spmatrix_t} )::spm_int_t
end

@cbindings libspm begin
    @cextern spmCheckAndCorrect( spm_in::Ptr{spmatrix_t}, spm_out::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmGenMat( type::spm_rhstype_t, nrhs::spm_int_t, spm::Ptr{spmatrix_t}, alpha::Ptr{Cvoid}, seed::Culonglong, A::Ptr{Cvoid}, lda::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmGenVec( type::spm_rhstype_t, spm::Ptr{spmatrix_t}, alpha::Ptr{Cvoid}, seed::Culonglong, x::Ptr{Cvoid}, incx::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmGenRHS( type::spm_rhstype_t, nrhs::spm_int_t, spm::Ptr{spmatrix_t}, opt_X::Ptr{Cvoid}, opt_ldx::spm_int_t, B::Ptr{Cvoid}, ldb::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmCheckAxb( eps::Cdouble, nrhs::spm_int_t, spm::Ptr{spmatrix_t}, opt_X0::Ptr{Cvoid}, opt_ldx0::spm_int_t, B::Ptr{Cvoid}, ldb::spm_int_t, X::Ptr{Cvoid}, ldx::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmExtractLocalRHS( nrhs::spm_int_t, spm::Ptr{spmatrix_t}, Bg::Ptr{Cvoid}, ldbg::spm_int_t, Bl::Ptr{Cvoid}, ldbl::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmReduceRHS( nrhs::spm_int_t, spm::Ptr{spmatrix_t}, Bg::Ptr{Cvoid}, ldbg::spm_int_t, Bl::Ptr{Cvoid}, ldbl::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmGatherRHS( nrhs::spm_int_t, spm::Ptr{spmatrix_t}, Bl::Ptr{Cvoid}, ldbl::spm_int_t, root::Cint, Bg::Ptr{Cvoid}, ldbg::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmIntConvert( n::spm_int_t, input::Ptr{Cint}, output::Ptr{spm_int_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmLoadDist( spm::Ptr{spmatrix_t}, filename::Cstring, comm::__get_mpi_type__() )::Cint
end

@cbindings libspm begin
    @cextern spmLoad( spm::Ptr{spmatrix_t}, filename::Cstring )::Cint
end

@cbindings libspm begin
    @cextern spmSave( spm::Ptr{spmatrix_t}, filename::Cstring )::Cint
end

@cbindings libspm begin
    @cextern spmReadDriver( driver::spm_driver_t, filename::Cstring, spm::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmReadDriverDist( driver::spm_driver_t, filename::Cstring, spm::Ptr{spmatrix_t}, comm::__get_mpi_type__() )::Cint
end

@cbindings libspm begin
    @cextern spmParseLaplacianInfo( filename::Cstring, flttype::Ptr{spm_coeftype_t}, dim1::Ptr{spm_int_t}, dim2::Ptr{spm_int_t}, dim3::Ptr{spm_int_t}, alpha::Ptr{Cdouble}, beta::Ptr{Cdouble}, dof::Ptr{spm_int_t} )::Cint
end

@cbindings libspm begin
    @cextern spm2Dense( spm::Ptr{spmatrix_t}, A::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmPrint( spm::Ptr{spmatrix_t}, f::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmPrintRHS( spm::Ptr{spmatrix_t}, nrhs::Cint, x::Ptr{Cvoid}, ldx::spm_int_t, stream::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmPrintInfo( spm::Ptr{spmatrix_t}, f::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmExpand( spm_in::Ptr{spmatrix_t}, spm_out::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmDofExtend( spm::Ptr{spmatrix_t}, type::Cint, dof::Cint, spm_out::Ptr{spmatrix_t} )::Cint
end

end #module
