#=

 @file spm.jl

 SPM julia wrapper

 @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Mathieu Faverge
 @author Lebdaoui selmane
 @date 2020-06-18

This file has been automatically generated with gen_wrappers.py

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

if spm_mpi_enabled
    using MPI
end

function __get_mpi_type__()
     if !spm_mpi_enabled
         return Cint
     elseif sizeof(MPI.MPI_Comm) == sizeof(Clong)
         return Clong
     elseif sizeof(MPI.MPI_Comm) == sizeof(Cint)
         return Cint
     end
     return Cvoid
end

@cstruct spmatrix_t {
    mtxtype::spm_mtxtype_t
    flttype::spm_coeftype_t
    fmttype::spm_fmttype_t
    gN::spm_int_t
    n::spm_int_t
    gnnz::spm_int_t
    nnz::spm_int_t
    gNexp::spm_int_t
    nexp::spm_int_t
    gnnzexp::spm_int_t
    nnzexp::spm_int_t
    dof::spm_int_t
    dofs::Ptr{spm_int_t}
    layout::spm_layout_t
    colptr::Ptr{spm_int_t}
    rowptr::Ptr{spm_int_t}
    loc2glob::Ptr{spm_int_t}
    values::Ptr{Cvoid}
    glob2loc::Ptr{spm_int_t}
    clustnum::Cint
    clustnbr::Cint
    comm::__get_mpi_type__()
}

@cbindings libspm begin
    @cextern spmInit( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmAlloc( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmExit( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmCopy( spm::Ptr{spmatrix_t} )::spmatrix_t
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
    @cextern spmInitDist( spm::Ptr{spmatrix_t}, comm::__get_mpi_type__() )::Cvoid
end

@cbindings libspm begin
    @cextern spmScatter( spm::Ptr{spmatrix_t}, n::spm_int_t, loc2glob::Ptr{spm_int_t}, distByColumn::Cint, root::Cint, comm::__get_mpi_type__() )::spmatrix_t
end

@cbindings libspm begin
    @cextern spmGather( spm::Ptr{spmatrix_t}, root::Cint )::spmatrix_t
end

@cbindings libspm begin
    @cextern spmNorm( ntype::spm_normtype_t, spm::Ptr{spmatrix_t} )::Cdouble
end

@cbindings libspm begin
    @cextern spmMatVec( trans::spm_trans_t, alpha::Cdouble, spm::Ptr{spmatrix_t}, x::Ptr{Cvoid}, beta::Cdouble, y::Ptr{Cvoid} )::Cint
end

@cbindings libspm begin
    @cextern spmMatMat( trans::spm_trans_t, n::spm_int_t, alpha::Cdouble, A::Ptr{spmatrix_t}, B::Ptr{Cvoid}, ldb::spm_int_t, beta::Cdouble, C::Ptr{Cvoid}, ldc::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmScalMatrix( alpha::Cdouble, spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmScalVector( flt::spm_coeftype_t, alpha::Cdouble, n::spm_int_t, x::Ptr{Cvoid}, incx::spm_int_t )::Cvoid
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
    @cextern spmGenRHS( type::spm_rhstype_t, nrhs::spm_int_t, spm::Ptr{spmatrix_t}, x::Ptr{Cvoid}, ldx::spm_int_t, b::Ptr{Cvoid}, ldb::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmCheckAxb( eps::Cdouble, nrhs::spm_int_t, spm::Ptr{spmatrix_t}, x0::Ptr{Cvoid}, ldx0::spm_int_t, b::Ptr{Cvoid}, ldb::spm_int_t, x::Ptr{Cvoid}, ldx::spm_int_t )::Cint
end

@cbindings libspm begin
    @cextern spmIntConvert( n::spm_int_t, input::Ptr{Cint} )::spm_int_t
end

@cbindings libspm begin
    @cextern spmLoad( spm::Ptr{spmatrix_t}, infile::Ptr{Cvoid} )::Cint
end

@cbindings libspm begin
    @cextern spmSave( spm::Ptr{spmatrix_t}, outfile::Ptr{Cvoid} )::Cint
end

@cbindings libspm begin
    @cextern spmReadDriver( driver::spm_driver_t, filename::Cstring, spm::Ptr{spmatrix_t} )::Cint
end

@cbindings libspm begin
    @cextern spmParseLaplacianInfo( filename::Cstring, flttype::Ptr{spm_coeftype_t}, dim1::Ptr{spm_int_t}, dim2::Ptr{spm_int_t}, dim3::Ptr{spm_int_t}, alpha::Ptr{Cdouble}, beta::Ptr{Cdouble}, dof::Ptr{spm_int_t} )::Cint
end

@cbindings libspm begin
    @cextern spm2Dense( spm::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmPrint( spm::Ptr{spmatrix_t}, f::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmPrintInfo( spm::Ptr{spmatrix_t}, f::Ptr{Cvoid} )::Cvoid
end

@cbindings libspm begin
    @cextern spmExpand( spm_in::Ptr{spmatrix_t}, spm_out::Ptr{spmatrix_t} )::Cvoid
end

@cbindings libspm begin
    @cextern spmDofExtend( spm::Ptr{spmatrix_t}, type::Cint, dof::Cint )::spmatrix_t
end

end   #module