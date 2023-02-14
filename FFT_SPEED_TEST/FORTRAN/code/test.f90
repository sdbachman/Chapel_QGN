program main

use ISO_C_BINDING
use :: MKL_DFTI
implicit none ! Applies to all contained subroutines

integer, parameter :: dp = C_DOUBLE !, dpc = C_DOUBLE_COMPLEX
real(dp), parameter :: pi = 4._dp*atan(1.0_dp)
REPL1
REPL2
REPL3

!integer, parameter :: nx = (3*128)/2
!integer, parameter :: ny = (3*128)/2
!integer, parameter :: nz = 8

integer :: i, k
real(dp) :: q0(nx,ny,nz) ! Initial condition
! Generic holders for FFT inputs
real(dp) :: grid(nx,ny,nz), gridE(nx*ny*nz)
Equivalence (grid, gridE)
! Generic holders for FFT outputs
complex(dp) :: spec(nx/2 +1,ny,nz), specE((nx/2 +1)*ny*nz)
Equivalence (spec, specE)

! Some MKL vars
integer :: Dims(2)
integer :: real_strides(3), cplx_strides(3)
integer :: zeroInd, maxInd ! index of zero and 1st baroclinic eigenvalues
integer :: Status
! MKL_DFTI plans
type(DFTI_DESCRIPTOR), pointer :: g2s_q, s2g_q, g2s_b, s2g_b

! Timing
integer :: start, finish, xx
real(dp) :: rate

! For writing the ICs
character(len=1024) :: filename


call RANDOM_NUMBER(q0)
do i=1,nz
    q0(:,:,i) = 1.E-6_dp*(q0(:,:,i) - sum(q0(:,:,i))/real(nx*ny,dp))
end do
do i=1,nx
    q0(i,:,:) = q0(i,:,:) + 1.E-5*sin(16._dp*pi*real(i-1,dp)/real(nx,dp))
end do

do k = 1,nz

if (k < 10) then
write (filename, "(A10,I1,A4)") "test_grid0", k, ".dat"
open(unit=10,file=trim(filename),access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*nx*ny)
    write(10,REC=1) q0(:,:,k)
    close(10)
else
write (filename, "(A9,I2,A4)") "test_grid", k, ".dat"
open(unit=10,file=trim(filename),access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=8*nx*ny)
    write(10,REC=1) q0(:,:,k)
    close(10)
end if



end do


! Set up MKL DFTI plans
    include "MKL_DFTI_Setup.f90"

! Transform the initial condition

call system_clock(start, rate)
do k = 1,10000
    grid = q0
    Status = DftiComputeForward(g2s_q, gridE, specE)
    !print *, Status
end do
call system_clock(finish)

write(*, *) (finish-start)/rate

end program main

