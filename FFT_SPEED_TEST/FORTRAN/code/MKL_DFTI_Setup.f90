Dims = [nx, ny]
real_strides(1) = 0
real_strides(2) = 1
real_strides(3) = nx
cplx_strides(1) = 0
cplx_strides(2) = 1
cplx_strides(3) = nx/2 + 1
!print *,'---------------------------------------------------------------'
!print *,' Printing a bunch of status outputs for MKL DFTI setup '
! Configure forward, real-to-complex fft, everything but b
Status = DftiCreateDescriptor(g2s_q, DFTI_DOUBLE, DFTI_REAL, 2, Dims)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_NUMBER_OF_TRANSFORMS, nz)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_INPUT_DISTANCE, nx*ny)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_OUTPUT_DISTANCE, (nx/2+1)*ny)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_INPUT_STRIDES, real_strides)
!print *, Status
Status = DftiSetValue(g2s_q, DFTI_OUTPUT_STRIDES, cplx_strides)
!print *, Status
Status = DftiCommitDescriptor(g2s_q)
!print *, Status
! Configure backward, complex-to-real fft, everything but b
Status = DftiCreateDescriptor(s2g_q, DFTI_DOUBLE, DFTI_REAL, 2, Dims)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_NUMBER_OF_TRANSFORMS, nz)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_INPUT_DISTANCE, (nx/2 +1)*ny)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_OUTPUT_DISTANCE, nx*ny)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_INPUT_STRIDES, cplx_strides)
!print *, Status
Status = DftiSetValue(s2g_q, DFTI_OUTPUT_STRIDES, real_strides)
!print *, Status
Status = DftiCommitDescriptor(s2g_q)
!print *, Status
!print *,'---------------------------------------------------------------'
