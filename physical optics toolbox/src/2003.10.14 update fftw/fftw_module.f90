! MODULE to make transparent calls to fftw, using the same interface as for Benson's
! fft calls
MODULE fftw_module

PUBLIC

    !! fftw constants copied from fftw3.f
      INTEGER, PARAMETER :: FFTW_R2HC = 0
      INTEGER, PARAMETER :: FFTW_HC2R = 1
      INTEGER, PARAMETER :: FFTW_DHT = 2
      INTEGER, PARAMETER :: FFTW_REDFT00=3
      INTEGER, PARAMETER :: FFTW_REDFT01=4
      INTEGER, PARAMETER :: FFTW_REDFT10=5
      INTEGER, PARAMETER :: FFTW_REDFT11=6
      INTEGER, PARAMETER :: FFTW_RODFT00=7
      INTEGER, PARAMETER :: FFTW_RODFT01=8
      INTEGER, PARAMETER :: FFTW_RODFT10=9
      INTEGER, PARAMETER :: FFTW_RODFT11=10
      INTEGER, PARAMETER :: FFTW_FORWARD = -1
      INTEGER, PARAMETER :: FFTW_BACKWARD=+1
      INTEGER, PARAMETER :: FFTW_MEASURE=0
      INTEGER, PARAMETER :: FFTW_DESTROY_INPUT=1
      INTEGER, PARAMETER :: FFTW_UNALIGNED=2
      INTEGER, PARAMETER :: FFTW_CONSERVE_MEMORY=4
      INTEGER, PARAMETER :: FFTW_EXHAUSTIVE=8
      INTEGER, PARAMETER :: FFTW_PRESERVE_INPUT=16
      INTEGER, PARAMETER :: FFTW_PATIENT=32
      INTEGER, PARAMETER :: FFTW_ESTIMATE=64
      INTEGER, PARAMETER :: FFTW_ESTIMATE_PATIENT=128
      INTEGER, PARAMETER :: FFTW_BELIEVE_PCOST=256
      INTEGER, PARAMETER :: FFTW_DFT_R2HC_ICKY=512
      INTEGER, PARAMETER :: FFTW_NONTHREADED_ICKY=1024
      INTEGER, PARAMETER :: FFTW_NO_BUFFERING=2048
      INTEGER, PARAMETER :: FFTW_NO_INDIRECT_OP=4096
      INTEGER, PARAMETER :: FFTW_ALLOW_LARGE_GENERIC=8192
      INTEGER, PARAMETER :: FFTW_NO_RANK_SPLITS=16384
      INTEGER, PARAMETER :: FFTW_NO_VRANK_SPLITS=32768
      INTEGER, PARAMETER :: FFTW_NO_VRECURSE=65536
      INTEGER, PARAMETER :: FFTW_NO_SIMD=131072


END MODULE fftw_module
