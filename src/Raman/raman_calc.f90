PROGRAM Raman_Calc

    USE map_data
    USE time_data
    USE freq_data
    USE hist_data
    USE input_module
    USE cli_data
    USE CLI
    IMPLICIT NONe
    INTEGER, PARAMETER :: nperchunk = 1000
    INTEGER :: chunk, iper, ioh, nchunks

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: w01
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mu01, eOH, a_para, a_perp
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: z0
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: efield
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: tcf, vv_tcf_tot, vh_tcf_tot

    CALL Read_CLI_Arguments

! *********************************************************************
! I.  Read The Input File, and set the data for the calculations.
!     Also ALLOCATE the relevant variables.
! *********************************************************************

    CALL Read_Input

    CALL Apply_CLI_Args


    nchunks = ceiling(real(noh)/real(nperchunk))

    ALLOCATE(w01(nperchunk,ntimes)); ALLOCATE(mu01(nperchunk,ntimes,3))
    ALLOCATE(eOH(nperchunk,ntimes,3))
    ALLOCATE(z0(nperchunk,3)); ALLOCATE(efield(ntimes))

    ALLOCATE(a_para(nperchunk, ntimes, 3)); ALLOCATE(a_perp(nperchunk, ntimes, 3))

    ALLOCATE(tcf(0:ncorr)); ALLOCATE(vv_tcf_tot(0:ncorr)); ALLOCATE(vh_tcf_tot(0:ncorr))
    vv_tcf_tot = DCMPLX(0.0,0.0); vh_tcf_tot = DCMPLX(0.0,0.0)
    w01_avg  = 0.0d0; w01_sq_avg = 0.0d0

! *********************************************************************
!   II. Calculate TCFS
! *********************************************************************

    DO chunk=1, nchunks
        WRITE(*,*) 'Calculating TCFs for chunk ', chunk
        w01 = 0.0; mu01 = 0.0; eOH = 0.0; a_para = 0.0; a_perp = 0.0
        DO iper=1, nperchunk
            ioh = (chunk-1)*nperchunk + iper
            if (ioh > noh) EXIT
            CALL Read_Field_File(ioh, efield(:), eOH(iper,:,:), z0(iper,:))
            WRITE(*,*) 'READ'
            CALL Get_w01(efield(:),  w01(iper,:))
            WRITE(*,*) 'READ1'
            CALL Get_01_Dipole(efield(:), w01(iper,:),  eOH(iper,:,:), mu01(iper,:,:))
            WRITE(*,*) 'READ2'
            CALL Get_Transition_Pol_Para_Perp(efield(:), w01(iper,:), eOH(iper,:,:), &
                                    & a_para(iper,:,:), a_perp(iper,:,:))
            WRITE(*,*) 'READ3'
        END DO

        WRITE(*,*) w01(1,1), mu01(1,1,1)

        DO ioh=(chunk-1)*nperchunk+1, min(chunk*nperchunk, noh)
                iper = ioh - (chunk-1)*nperchunk
                ! VV TCF
                CALL Calc_TCF(w01(iper,:), a_para(iper,:,:), tcf(:))
                vv_tcf_tot = vv_tcf_tot + tcf(:)
                ! VH TCF
                CALL Calc_TCF(w01(iper,:), a_perp(iper,:,:), tcf(:))
                vh_tcf_tot = vh_tcf_tot + tcf(:)
        ENDDO 
    END DO

    ! Normalize by the number of oh groups
    vv_tcf_tot = vv_tcf_tot/DCMPLX(DFLOAT(noh),0d0)
    vh_tcf_tot = vh_tcf_tot/DCMPLX(DFLOAT(noh),0d0)
    w01_avg = w01_avg/DFLOAT(noh*ntimes)
    w01_sq_avg = w01_sq_avg/DFLOAT(noh*ntimes)

! *********************************************************************
! III. Calculate the Raman Spectra
! *********************************************************************

    CALL Spec_Calc(vv_tcf_tot, vh_tcf_tot)

! *********************************************************************
! IV. Cleanup Calculation
! *********************************************************************

DEALLOCATE(z0, efield)
DEALLOCATE(w01, mu01, eOH,  a_para, a_perp, tcf, vv_tcf_tot, vh_tcf_tot)

END PROGRAM Raman_Calc