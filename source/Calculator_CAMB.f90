    !Use CAMB
    module Calculator_CAMB
    use CosmologyTypes
    use CosmoTheory
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
        outNone,  CAMBdata, NonLinear_Pk, Nonlinear_lens, CAMB_TransfersToPowers, initial_nummodes, &
        initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
        highL_unlensed_cl_template, nthermo_derived, TInitialPowerLaw, &
        threadnum, version, tensor_param_rpivot
    use NonLinear, only : halofit_default, THalofit
    use config
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
    use, intrinsic :: ieee_arithmetic
    implicit none
    private

    Type, extends(TTheoryIntermediateCache) :: CAMBTransferCache
        Type (CAMBdata) :: State
    end Type CAMBTransferCache

    Type, extends(TCosmologyCalculator) :: CAMB_Calculator
        integer :: ncalls = 0
        integer :: nerrors = 0
        logical :: CAMB_timing = .false.
        real(mcp) :: k_eta_max_scalar = -1._mcp
        logical :: accurate_BB =.false.
        integer :: halofit_version = halofit_default
        character(LEN=:), allocatable :: dark_energy_model
        type(CAMBParams)  CAMBP
        character(LEN=:), allocatable :: highL_theory_cl_template_file
        real(mcp), allocatable :: highL_lensedCL_template(:,:)
        Type(CAMBTransferCache) :: DefaultInstance
        class(CAMBdata), pointer :: CurrentState
    contains
    !New
    procedure :: SetCurrentPoint => CAMBcalc_SetCurrentPoint
    procedure :: CMBToCAMB => CAMBCalc_CMBToCAMB
    procedure :: SetDerived => CAMBCalc_SetDerived
    procedure :: SetPowersFromCAMB => CAMBCalc_SetPowersFromCAMB
    procedure :: InitCAMB => CAMBCalc_InitCAMB
    procedure :: InitCAMBParams => CAMBCalc_InitCAMBParams
    procedure :: SetCAMBInitPower => CAMBCalc_SetCAMBInitPower
    procedure :: SetPkFromCAMB => CAMBCalc_SetPkFromCAMB
    procedure :: GetNLandRatios => CAMBCalc_GetNLandRatios
    !Overridden inherited
    procedure :: ReadParams => CAMBCalc_ReadParams
    procedure :: InitForLikelihoods => CAMBCalc_InitForLikelihoods
    procedure :: BAO_D_v => CAMBCalc_BAO_D_v
    procedure :: AngularDiameterDistance => CAMBCalc_AngularDiameterDistance
    procedure :: ComovingRadialDistance => CAMBCalc_ComovingRadialDistance
    procedure :: ComovingRadialDistanceArr => CAMBCalc_ComovingRadialDistanceArr
    procedure :: AngularDiameterDistance2 => CAMBCalc_AngularDiameterDistance2
    procedure :: LuminosityDistance => CAMBCalc_LuminosityDistance
    procedure :: Hofz => CAMBCalc_Hofz
    procedure :: CMBToTheta => CAMBCalc_CMBToTheta
    procedure :: GetNewPowerData => CAMBCalc_GetNewPowerData
    procedure :: GetNewTransferData => CAMBCalc_GetNewTransferData
    procedure :: GetCosmoTheoryForImportance => CAMBCalc_GetTheoryForImportance
    procedure :: GetZreFromTau  => CAMBCalc_GetZreFromTau
    procedure :: GetOpticalDepth => CAMBCalc_GetOpticalDepth
    procedure :: SetBackgroundTheoryData => CAMBCalc_SetBackgroundTheoryData
    procedure :: SetParamsForBackground => CAMBCalc_SetParamsForBackground
    procedure :: VersionTraceOutput => CAMBCalc_VersionTraceOutput
    procedure, private :: LoadFiducialHighLTemplate
    end type CAMB_Calculator


    !    integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
    !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT


    public CAMB_Calculator
    contains

    subroutine CAMBCalc_CMBToCAMB(this,CMB,P)
    use DarkEnergyInterface
    use Reionization
    use model
    use constants, only : dl, default_nnu,delta_mnu21,delta_mnu31,mnu_min_normal
    use lensing, only : ALens_Fiducial
    use MassiveNu, only : sum_mnu_for_m1
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    type(CAMBParams)  P
    real(dl) neff_massive_standard, mnu, m1, m3, normal_frac
    real(dl), external :: Newton_raphson

    P = this%CAMBP
    P%ombh2 = CMB%ombh2
    P%omnuh2 = CMB%omnuh2
    P%omch2 = CMB%omch2
    P%omk = CMB%omk
    P%H0 = CMB%H0
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
        Reion%redshift= CMB%zre
        Reion%delta_redshift = CMB%zre_delta
    class default
        call MpiStop('Unknown reionizxation model')
    end select
    select type(DE => P%DarkEnergy)
    class is (TDarkEnergyEqnOfState)
        DE%w_lam = CMB%w
        DE%wa = CMB%wa
    class default
        call MpiStop('Unknown dark energy model')
    end select
    P%ALens = CMB%ALens
    ALens_Fiducial = CMB%ALensf
    P%InitialConditionVector(initial_iso_CDM) = &
        sign(sqrt(abs(CMB%iso_cdm_correlated) /(1-abs(CMB%iso_cdm_correlated))),CMB%iso_cdm_correlated)
    P%Num_Nu_Massive = 0
    P%Nu_mass_numbers = 0
    P%Num_Nu_Massless = CMB%nnu
    P%share_delta_neff = .false.

    if (CMB%omnuh2>0) then
        call P%SetNeutrinoHierarchy(CMB%omnuh2, CMB%omnuh2_sterile, CMB%nnu, &
            CosmoSettings%neutrino_hierarchy, CosmoSettings%num_massive_neutrinos)
    end if

    P%YHe = CMB%YHe
#ifdef COSMOREC
    if (CMB%fdm/=0._mcp) P%Recomb%runmode = 3
    P%Recomb%fdm = CMB%fdm * 1e-23_mcp
#else
    if (CMB%fdm/=0._mcp) call MpiStop('Compile with CosmoRec to use fdm')
#endif
    call this%SetCAMBInitPower(P,CMB)

    end subroutine CAMBCalc_CMBToCAMB

    subroutine CAMBCalc_SetDerived(this, Info, Theory)
    class(CAMB_Calculator) :: this
    class(TTheoryIntermediateCache) :: Info
    Class(TCosmoTheoryPredictions) Theory
    integer noutputs, i

    select type (Info)
    class is (CAMBTransferCache)
        noutputs = size(Info%State%CP%z_outputs)
        Theory%numderived = nthermo_derived + noutputs*2
        if (Theory%numderived > max_derived_parameters) &
            call MpiStop('numderived > max_derived_parameters: increase in CosmologyTypes.f90')
        Theory%derived_parameters(1:nthermo_derived) = Info%State%ThermoDerivedParams(1:nthermo_derived)
        do i=1, noutputs
            !Theory%derived_parameters(nthermo_derived+(i-1)*3+1) = BackgroundOutputs%rs_by_D_v(i)
            !now use Hubble paramter in normal units and DM, comoving angular diameter distance
            Theory%derived_parameters(nthermo_derived+(i-1)*2+1) = Info%State%BackgroundOutputs%H(i)*const_c/1e3_mcp
            Theory%derived_parameters(nthermo_derived+(i-1)*2+2) = &
                Info%State%BackgroundOutputs%DA(i)*(1+Info%State%CP%z_outputs(i))
            !Theory%derived_parameters(nthermo_derived+(i-1)*4+4) = (1+BackgroundOutputs%z_outputs(i))* &
            !    BackgroundOutputs%DA(i) * BackgroundOutputs%H(i) !F_AP parameter
        end do
    end select

    end subroutine CAMBCalc_SetDerived

    subroutine CAMBCalc_SetParamsForBackground(this,CMB, Info)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Type(CAMBParams)  P
    class(TTheoryIntermediateCache), target :: Info

    select type (Info)
    class is (CAMBTransferCache)
        !set background parameters, but don't calculate thermal history
        call this%CMBToCAMB(CMB, P)
        call Info%State%SetParams(P)
        this%CurrentState => Info%State
    end select

    end subroutine CAMBCalc_SetParamsForBackground

    subroutine CAMBCalc_SetBackgroundTheoryData(this, CMB, Info, Theory,error)
    use cambmain, only: initvars
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TTheoryIntermediateCache) :: Info
    Class(TCosmoTheoryPredictions) Theory
    integer error

    select type (Info)
    class is (CAMBTransferCache)
        call InitVars(Info%State) !calculate thermal history, e.g. z_drag etc.
        if (global_error_flag/=0) then
            error=global_error_flag
            return
        end if
        call this%SetDerived(Info, Theory)
    end select

    end subroutine CAMBCalc_SetBackgroundTheoryData

    subroutine CAMBCalc_GetNewTransferData(this, CMB,Info,Theory,error)
    use InitialPower
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error
    type(CAMBParams)  P
    real(mcp) time

    allocate(Info, source=this%DefaultInstance)
    select type (Info)
    class is (CAMBTransferCache)

        call this%CMBToCAMB(CMB, P)
        this%CurrentState => Info%State

        if (Feedback > 1) write (*,*) 'Calling CAMB'
        Threadnum =num_threads

        if (this%CAMB_timing) time = TimerTime()
        call CAMB_GetResults(Info%State, P, error, onlytransfer=.true.)
        if (this%CAMB_timing) call Timer('GetTransfers', time)
    class default
        call MpiStop('CAMB_GetNewTransferData with wrong TTheoryIntermediateCache type')
    end select
    if (error==0) then
        call this%SetDerived(Info,Theory)
    else
        if (stop_on_error) call MpiStop('CAMB error '//trim(global_error_message))
        if (Feedback > 0) write(*,*) 'CAMB returned error '//trim(global_error_message)
        this%nerrors=this%nerrors+1
    end if
    this%ncalls=this%ncalls+1
    if (mod(this%ncalls,100)==0 .and. LogFile%Opened()) then
        write(LogFile%unit,'("CAMB called ",I0," times; ",I0," errors")') this%ncalls, this%nerrors
    end if
    if (Feedback > 1) write (*,*) 'CAMB done'

    end subroutine CAMBCalc_GetNewTransferData

    subroutine CAMBCalc_GetNewPowerData(this, CMB, Info, Theory, error)
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer :: error,i,j

    select type (Info)
    class is (CAMBTransferCache)
        call this%SetCAMBInitPower(Info%State%CP,CMB)
        call CAMB_TransfersToPowers(Info%State)
        !this sets slow CAMB params correctly from value stored in Transfers
        if (global_error_flag/=0) then
            error=global_error_flag
            return
        end if
        !JD 08/13 added so we dont have to fill Cls unless using CMB
        if(CosmoSettings%use_CMB)then
            call this%SetPowersFromCAMB(CMB,Info%State, Theory)
            do i=1, min(3,CosmoSettings%num_cls)
                if (CosmoSettings%cl_lmax(i,i)>0) then
                    if (any(Theory%cls(i,i)%Cl(:) < 0 )) then
                        error = 1
                        !write(*,*) '1) FIELD IS', i !, Theory%cls(i,i)%Cl(:)!!NR
                        call MpiStop('Calculator_CAMB: negative C_l (could edit to silent error here)')
                        return
                    end if
                end if
                do j= i, 1, -1
                    if (CosmoSettings%cl_lmax(i,j)>0) then
                        if ( any(ieee_is_nan(Theory%cls(i,j)%Cl))) then
                            error=1
                            write(*,*) 'WARNING: NaN CL?', i, j
                            return
                        end if
                    end if
                end do
            end do
        end if

        !redshifts are in increasing order, so last index is redshift zero
        if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
            Theory%sigma_8 = Info%State%MT%sigma_8(size(Info%State%MT%sigma_8))
        else
            Theory%sigma_8 = 0
        end if

        if (CosmoSettings%Use_LSS) then
            call this%SetPkFromCAMB(Info%State,Theory,error)
            if (error/=0) then
                write(*,*) 'WARNING: NaN PK?'
                return
            end if
        end if
    end select

    end subroutine CAMBCalc_GetNewPowerData

    subroutine CAMBCalc_GetTheoryForImportance(this, CMB, Theory, error)
    use config, only : ThreadNum
    class(CAMB_Calculator), target :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    class(CAMBTransferCache), pointer :: Info
    integer error,i,j
    logical :: DoCls, DoPk
    type(CAMBParams)  P

    error = 0
    DoCls = this%ImportanceOptions%redo_cls
    DoPk = this%ImportanceOptions%redo_pk
    Info => this%DefaultInstance

    if (DoCls .or. DoPk) then
        call this%CMBToCAMB(CMB, P)
        Threadnum =num_threads

        P%WantCls = DoCls

        if (.not. DoPk .and. .not. (CosmoSettings%CMB_Lensing .and. &
            CosmoSettings%use_nonlinear_lensing)) P%WantTransfer = .false.

        if (this%CAMB_timing) call Timer()
        if (Feedback > 1) write (*,*) 'Calling CAMB'
        call CAMB_GetResults(Info%State, P, error, onlytransfer=.true.)

        if (Feedback > 1) write (*,*) 'CAMB Done'
        if (this%CAMB_timing) call Timer('CAMB_GetTransfers')

        if (error==0) then
            call this%SetCAMBInitPower(Info%State%CP,CMB)
            call CAMB_TransfersToPowers(Info%State)
            error=global_error_flag
        end if
    else
        call this%GetNewBackgroundData(CMB,Info,Theory,error)
    end if

    if (DoCls .and. error==0) then
        call this%SetPowersFromCAMB(CMB, Info%State, Theory)
        if (any(Theory%cls(1,1)%Cl(:) < 0 )) then
            error = 1
            !write(*,*) '2) FIELD IS', i !!NR
            call MpiStop('Calculator_CAMB: negative C_l (could edit to silent error here)')
        end if
        do i=1, min(3,CosmoSettings%num_cls)
            if(error/=0) exit
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(ieee_is_nan(Theory%cls(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?'
                        exit
                    end if
                end if
            end do
        end do
    end if

    if (DoPK .and. error==0) then
        Theory%sigma_8 = Info%State%MT%sigma_8(size(Info%State%MT%sigma_8))
        call this%SetPkFromCAMB(Info%State,Theory,error)
        if (error/=0) write(*,*) 'WARNING: NaN PK?'
    end if

    if (error==0) call this%SetDerived(Info, Theory)

    call Info%Clear()

    end subroutine CAMBCalc_GetTheoryForImportance

    subroutine CAMBCalc_SetPowersFromCAMB(this,CMB, State, Theory)
    use constants
    use InitialPower
    use config,only : C_Temp, C_E, CT_Cross, CT_B, C_Cross, C_Phi, C_PhiTemp
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(CAMBdata) :: State
    class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), parameter :: cons =  (COBE_CMBTemp*1e6)**2
    integer l
    real(mcp) :: highL_norm = 0
    real(mcp) lens_recon_scale, rms
    integer :: i,j, lmx, lmaxCL, C_V = 4 !!NR
    integer, save, allocatable :: indicesS(:,:), indicesT(:,:)
    real(mcp), dimension(:), allocatable :: newCl_E, newCl_TE, newCl_B, CLBB, Cl_VV !!NR
    integer :: i_F, index_fw3j, myunit
    real(mcp) :: prefac_cl2dl, currclEE, currclTE, currclBB, currclEE_p1, currclBB_p1, currclEE_m1, currclBB_m1, currclBB_p2, currclBB_m2 !!NR
    !logical :: CLfile = .true. !!NR flag to print spectra on file cl_out.txt

    if (.not. allocated(indicesS)) then
        allocate(indicesS(4,4))   !!NR
        allocate(indicesT(4,4))   !!NR
        indicesS = 0
        indicesS(1,1) = C_Temp
        indicesS(2,2) = C_E
        indicesS(4,4) = C_V  !!NR
        indicesT = indicesS
        indicesT(2,1) = CT_Cross
        indicesT(3,3) = CT_B
        indicesS(2,1) = C_Cross
    end if

    lens_recon_scale = CMB%InitPower(Aphiphi_index)
    
    !!NR
    if(CosmoSettings%fw3jread .eq. .false.) then
        if (allocated(CosmoSettings%W1)) deallocate(CosmoSettings%W1)
        allocate(CosmoSettings%W1(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W1m1)) deallocate(CosmoSettings%W1m1)
        allocate(CosmoSettings%W1m1(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W1p1)) deallocate(CosmoSettings%W1p1)
        allocate(CosmoSettings%W1p1(CosmoSettings%lmax_dc))

        if (allocated(CosmoSettings%W2)) deallocate(CosmoSettings%W2)   !!NR aggiunto i Wigner anche per costruire VV
        allocate(CosmoSettings%W2(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W2p1)) deallocate(CosmoSettings%W2p1)
        allocate(CosmoSettings%W2p1(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W2m1)) deallocate(CosmoSettings%W2m1)
        allocate(CosmoSettings%W2m1(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W2p2)) deallocate(CosmoSettings%W2p2)
        allocate(CosmoSettings%W2p2(CosmoSettings%lmax_dc))
        if (allocated(CosmoSettings%W2m2)) deallocate(CosmoSettings%W2m2)
        allocate(CosmoSettings%W2m2(CosmoSettings%lmax_dc))
        
        open (newunit=myunit, file=trim(CosmoSettings%infile_fw3j_ee_bb), status='old', action='read')
        do i_F=1,CosmoSettings%lmax_dc,1
            read(myunit,*) CosmoSettings%W1(i_F), CosmoSettings%W1p1(i_F), CosmoSettings%W1m1(i_F)
        enddo
        close(myunit)

        open (newunit=myunit, file=trim(CosmoSettings%infile_fw3j_vv), status='old', action='read')   !!NR anche fw3j_vv
        do i_F=1,CosmoSettings%lmax_dc,1
            read(myunit,*) CosmoSettings%W2(i_F), CosmoSettings%W2p1(i_F), CosmoSettings%W2m1(i_F), CosmoSettings%W2p2(i_F), CosmoSettings%W2m2(i_F)
        enddo
        close(myunit)
        CosmoSettings%fw3jread = .true.
    endif
    !!NR

    do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%Cls(i,j)%CL)
                    if (indicesT(i,j)==0) then
                        CL=0
                    else
                        if (CosmoSettings%CMB_Lensing) then
                            CL(2:lmx) = cons*State%CLData%Cl_lensed(2:lmx, indicesT(i,j))
                        else
                            if (indicesS(i,j)/=0) then
                                CL(2:lmx) = cons*State%CLData%Cl_Scalar(2:lmx, indicesS(i,j))
                            else
                                CL=0
                            end if
                        end if
                        if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                            if (highL_norm ==0) & !normally normalize off TT
                                & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                            CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                        end if
                        if (CosmoSettings%compute_tensors) then
                            lmx = min(lmx,CosmoSettings%lmax_tensor)
                            CL(2:lmx) =  CL(2:lmx) + cons*State%CLData%Cl_tensor(2:lmx, indicesT(i,j))
                        end if
                    end if
                end associate
            end if
        end do
    end do

    if (CosmoSettings%use_lensing_potential) then
        !CMB lensing potential
        !in camb Cphi is l^4 C_l, we want [l(l+1)]^2Cphi/2pi
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_Phi))
        if (lmx/=0) then
            if (.not. CosmoSettings%CMB_lensing) call MpiStop('Must have lensing on to use lensing potential')
            associate(CL=> Theory%Cls(CL_Phi,CL_Phi)%CL)
                do l=2, lmx
                    CL(L) =  State%CLData%Cl_scalar(L,C_Phi)*(real(l+1)**2/l**2)/twopi * lens_recon_scale
                end do
                CL(lmx+1:)=0
            end associate
        end if
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_T))
        if (lmx/=0) then
            !lensing-temp
            do l=2, lmx
                Theory%Cls(CL_phi,CL_T)%CL = State%CLData%Cl_scalar(l,C_PhiTemp)/real(l)**3 * sqrt(lens_recon_scale)
            end do
        end if
    end if

    if (CosmoSettings%CMB_Lensing .and. this%CAMBP%max_l>=2000) then
        !Get RMS deflection angle in arcmin
        rms=0
        do L=2, 2000
            rms = rms +  State%CLData%Cl_scalar(L,C_Phi)*(real(l+1)**2/l**2)/twopi*(L+0.5_mcp)/(L*(L+1))
        end do
        Theory%Lensing_rms_deflect = sqrt(rms)*180/pi*60
    else
        Theory%Lensing_rms_deflect = 0
    end if

    if (CosmoSettings%compute_tensors) then
        Theory%tensor_ratio_02 = State%CP%InitPower%TensorPower(0.002d0)/State%CP%InitPower%ScalarPower(0.002d0)
        Theory%tensor_AT = State%CP%InitPower%TensorPower(CosmoSettings%tensor_pivot_k)
        Theory%tensor_ratio_BB = State%CP%InitPower%TensorPower(0.01d0)/State%CP%InitPower%ScalarPower(0.01d0)
        Theory%tensor_ratio_C10 = State%CLData%Cl_tensor(10, 1)/State%CLData%Cl_scalar(10, 1)
    else
        Theory%tensor_ratio_02 = 0
        Theory%tensor_ratio_BB = 0
        Theory%tensor_ratio_C10 = 0
        Theory%tensor_AT = 0
    end if


!!NR - calcolo i cl modificati
    if (CosmoSettings%fw3jread .eq. .true.) then  
        ! if ( (abs(CMB%fbv) .lt. CosmoSettings%tol_dc) ) then 
        !     CMB%fbv = 0.0
        ! endif
        allocate(newCl_E(2:CosmoSettings%lmax_dc))
        newCl_E = 0

        allocate(newCl_TE(2:CosmoSettings%lmax_dc))
        newCl_TE = 0

        allocate(newCl_B(2:CosmoSettings%lmax_dc))
        newCl_B = 0

        allocate(CLBB(0:CosmoSettings%lmax_dc+2))
        CLBB = 0

        allocate(Cl_VV(2:CosmoSettings%lmax_dc))   !!NR back aggiunto campo VV
        Cl_VV = 0

        if (CosmoSettings%lmax_tensor .le. CosmoSettings%lmax_dc+2) then
            CLBB(2:CosmoSettings%lmax_tensor) = cons*State%CLData%Cl_tensor(2:CosmoSettings%lmax_tensor,CT_B)
        else
            CLBB(2:CosmoSettings%lmax_dc+2) = cons*State%CLData%Cl_tensor(2:CosmoSettings%lmax_dc+2,CT_B)
        endif

            do l = 2, CosmoSettings%lmax_dc, 1
                index_fw3j = l-1
                prefac_cl2dl =l*(l+1.0)/twopi
                
                currclEE = Theory%Cls(2,2)%CL(l)*twopi/l/(l+1.0)
                currclBB = CLBB(l)*twopi/l/(l+1.0)
                currclTE = Theory%Cls(2,1)%CL(l)*twopi/l/(l+1.0)

                currclEE_p1 = Theory%Cls(2,2)%CL(l+1)*twopi/((l+2.0)*(l+1.0))
                currclBB_p1 = CLBB(l+1)*twopi/((l+2.0)*(l+1.0))

                currclEE_m1 = Theory%Cls(2,2)%CL(l-1)*twopi/(l*(l-1.0))
                currclBB_m1 = CLBB(l-1)*twopi/(l*(l-1.0))

                currclBB_p2 = CLBB(l+2)*twopi/((l+2.0)*(l+3.0))   !!NR
                if (l .le. 2) then
                        currclBB_m2 = 0.0
                else
                        currclBB_m2 = CLBB(l-2)*twopi/((l-2.0)*(l-1.0))   !!NR
                endif
              
                newCl_E(l) = (currclEE - 1/(2*twopi) * (CMB%fbv+CMB%fbe) * currclEE + 1/(2*twopi) * CMB%fbv * &
         (CosmoSettings%W1(index_fw3j)*currclEE+CosmoSettings%W1p1(index_fw3j)*currclBB_p1+CosmoSettings%W1m1(index_fw3j)*currclBB_m1)) *prefac_cl2dl

                newCl_B(l) = (currclBB - 1/(2*twopi) * (CMB%fbv+CMB%fbe) * currclBB + 1/(2*twopi) * CMB%fbv * & 
         (CosmoSettings%W1(index_fw3j)*currclBB+CosmoSettings%W1p1(index_fw3j)*currclEE_p1+CosmoSettings%W1m1(index_fw3j)*currclEE_m1)) *prefac_cl2dl

                newCl_TE(l) = (currclTE - 0.5/(2*twopi) * (CMB%fbv+CMB%fbe) * currclTE) *prefac_cl2dl
                
                Cl_VV(l) = CMB%fbe/(twopi/2) * (CosmoSettings%W2p2(index_fw3j)*currclBB_p2 + CosmoSettings%W2p1(index_fw3j)*currclEE_p1 + &
         CosmoSettings%W2(index_fw3j)*currclBB + CosmoSettings%W2m1(index_fw3j)*currclEE_m1 + CosmoSettings%W2m2(index_fw3j)*currclBB_m2) *prefac_cl2dl
                !!NR aggiunto spettro in VV
            enddo !closing l loop

            Theory%Cls(2,2)%CL(2:CosmoSettings%lmax_dc) = newCl_E 
            Theory%Cls(2,1)%CL(2:CosmoSettings%lmax_dc) = newCl_TE
            Theory%Cls(4,4)%CL(2:CosmoSettings%lmax_dc) = Cl_VV   !!NR
            if (allocated(Theory%Cls(3,3)%CL)) Theory%Cls(3,3)%CL(2:CosmoSettings%lmax_dc) = newCl_B + &
                CosmoSettings%alpha_dl*cons*State%CLData%Cl_lensed(2:CosmoSettings%lmax_dc, CT_B) !!NR alpha_dl is the delensing
             
            !write(*,*) 'fbe, fbv ', CMB%fbe, CMB%fbv
            if (CosmoSettings%CLfile) then !!NR
                open(107, file='cl_the.txt', status='replace')
                write(107,'("#",1A'//Trim(IntToStr(4))//'," ",*(A15))') 'L','TT','EE','TE','BB','VV'
                write(*,*) 'writing cl_the.txt THEORY' !!NR
                do i=2,5000
                write(107,'(1I6,5E20.12)') i,Theory%Cls(1,1)%CL(i), Theory%Cls(2,2)%CL(i), Theory%Cls(2,1)%CL(i), Theory%Cls(3,3)%CL(i), Theory%Cls(4,4)%CL(i)
                end do
                close(107)           
                CosmoSettings%CLfile = .false.
            endif
            deallocate(newCl_E,newCl_TE,newCl_B,CLBB,Cl_VV) !!NR

        endif
!!NR

    end subroutine CAMBCalc_SetPowersFromCAMB

    subroutine CAMBCalc_SetPkFromCAMB(this,State,Theory,error)
    use Transfer
    use camb, only : CP
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    class(CAMBdata) :: State
    integer :: error
    real(mcp), allocatable :: k(:), z(:), PK(:,:)
    integer zix,nz,nk, nR
    real(mcp), allocatable :: NL_Ratios(:,:)
    real(mcp) :: dR, R, minR
    integer i

    !Free theory arrays as they may resize between samples
    call Theory%FreePK()

    nz=CP%Transfer%PK_num_redshifts

    if (.not. allocated(Theory%growth_z)) allocate(Theory%growth_z, Theory%sigma8_z)
    call Theory%growth_z%InitForSize(nz)
    call Theory%sigma8_z%InitForSize(nz)
    allocate(z(nz))

    do zix=1,nz
        z(zix) = CP%Transfer%PK_redshifts(nz-zix+1)
        Theory%sigma8_z%F(zix) = State%MT%sigma_8(nz-zix+1)
        Theory%growth_z%F(zix) = State%MT%sigma2_vdelta_8(nz-zix+1)/State%MT%sigma_8(nz-zix+1)
    end do
    Theory%sigma8_z%X=z
    Theory%growth_z%X=z

    if (CosmoSettings%use_matterpower) then
        nk=State%MT%num_q_trans
        nz=CP%Transfer%PK_num_redshifts
        allocate(PK(nk,nz))
        allocate(k(nk))

        k = log(State%MT%TransferData(Transfer_kh,:,1))

        call Transfer_GetUnsplinedPower(State, State%MT, PK,transfer_power_var,transfer_power_var)
        PK = Log(PK)
        if (any(ieee_is_nan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK)
        call Theory%MPK%InitExtrap(k,z,PK, CosmoSettings%extrap_kmax)
    end if


    if (CosmoSettings%use_Weylpower) then
        call Transfer_GetUnsplinedPower(State,State%MT, PK,transfer_Weyl,transfer_Weyl,hubble_units=.false.)
        PK = Log(PK)
        if (any(ieee_is_nan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK_WEYL)
        call Theory%MPK_WEYL%InitExtrap(k,z,PK,CosmoSettings%extrap_kmax)
        ! Weyl density cross correlation:
        call Transfer_GetUnsplinedPower(State,State%MT, PK,transfer_Weyl,transfer_power_var,hubble_units=.false.)
        allocate(Theory%MPK_WEYL_CROSS)
        Theory%MPK_WEYL_CROSS%islog = .False.
        call Theory%MPK_WEYL_CROSS%InitExtrap(k,z,PK,CosmoSettings%extrap_kmax)
    end if

    if (CosmoSettings%use_SigmaR) then
        !Note R is in k*h units
        dR = log(1.2)/AccuracyLevel
        minR = 1/CP%Transfer%kmax
        nR = nint(log(150/minR)/dR) +1
        if (.not. allocated(Theory%Sigma_R)) allocate(Theory%Sigma_R)
        call Theory%Sigma_R%InitForSize(nR)
        do i=1, nR
            Theory%Sigma_R%X(i) = exp((i-1)*dR)*minR
        end do
        call Transfer_GetSigmaRArray(State, State%MT,Theory%Sigma_R%X, Theory%Sigma_R%F, &
            var1 = transfer_nonu,var2=transfer_nonu)
    end if

    if(CosmoSettings%use_nonlinear)then
        call this%GetNLandRatios(State,Theory,NL_Ratios,error)
        if(error/=0) return
    end if

    end subroutine CAMBCalc_SetPkFromCAMB


    subroutine CAMBCalc_GetNLandRatios(this,State,Theory,Ratios,error)
    use Transfer
    class(CAMB_Calculator) :: this
    class(CAMBdata) :: State
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    real(mcp), allocatable, intent(out) :: Ratios(:,:)
    Type(MatterPowerData) :: CPK
    real(mcp), allocatable :: PK(:,:)
    integer error,zix,nz

    CPK%num_k = Theory%MPK%nx
    CPK%num_z = Theory%MPK%ny

    !Allocate Theory arrays
    allocate(Theory%NL_MPK)
    allocate(Ratios(CPK%num_k,CPK%num_z))

    !fill PK with Linear MPK
    allocate(PK(CPK%num_k,CPK%num_z))
    PK=Theory%MPK%z

    allocate(CPK%matpower(CPK%num_k,CPK%num_z))
    allocate(CPK%ddmat(CPK%num_k,CPK%num_z))
    allocate(CPK%nonlin_ratio(CPK%num_k,CPK%num_z))
    allocate(CPK%log_kh(CPK%num_k))
    allocate(CPK%redshifts(CPK%num_z))
    CPK%log_kh = Theory%MPK%x
    CPK%redshifts = Theory%MPK%y
    CPK%matpower = PK

    !need splines to get nonlinear ratios
    call MatterPowerdata_getsplines(CPK)
    call this%CAMBP%NonLinearModel%GetNonLinRatios(State, CPK)
    Ratios = CPK%nonlin_ratio
    call MatterPowerdata_Free(CPK)

    PK = PK+2*log(Ratios)
    if (any(ieee_is_nan(PK))) then
        error = 1
        return
    end if
    call Theory%NL_MPK%InitExtrap(Theory%MPK%x,Theory%MPK%y,PK,CosmoSettings%extrap_kmax)

    if (allocated(Theory%MPK_WEYL)) then
        !Assume Weyl scales the same way under non-linear correction
        allocate(Theory%NL_MPK_WEYL)
        PK = Theory%MPK_WEYL%z + 2*log(Ratios)
        call Theory%NL_MPK_WEYL%InitExtrap(Theory%MPK%x,Theory%MPK%y,PK,CosmoSettings%extrap_kmax)
    end if

    if (allocated(Theory%MPK_WEYL_CROSS)) then
        !Assume Weyl scales the same way under non-linear correction
        allocate(Theory%NL_MPK_WEYL_CROSS)
        PK = Theory%MPK_WEYL_CROSS%z*Ratios**2
        Theory%NL_MPK_WEYL_CROSS%islog = .False.
        call Theory%NL_MPK_WEYL_CROSS%InitExtrap(Theory%MPK%x,Theory%MPK%y,PK,CosmoSettings%extrap_kmax)
    end if

    end subroutine CAMBCalc_GetNLandRatios

    subroutine CAMBCalc_InitCAMB(this,CMB,error, DoReion)
    class(CAMB_Calculator) :: this
    class(CMBParams), intent(in) :: CMB
    logical, optional, intent(in) :: DoReion
    logical WantReion
    type(CAMBParams)  P
    integer error

    if (present(DoReion)) then
        WantReion = DoReion
    else
        WantReion = .true.   !NR this has to be true
    end if

    call this%CMBToCAMB(CMB, P)
    call this%CurrentState%SetParams(P,error,WantReion)

    end subroutine CAMBCalc_InitCAMB

    function CAMBCalc_GetOpticalDepth(this,CMB) result(GetOpticalDepth)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) GetOpticalDepth
    type(CAMBParams)  P
    integer error
    Type(CAMBdata) State

    call this%CMBToCAMB(CMB, P)
    call State%SetParams(P,error)

    if (error/= 0) then
        GetOpticalDepth = -1
    else
        GetOpticalDepth = State%GetReionizationOptDepth()
    end if
    end function CAMBCalc_GetOpticalDepth

    function CAMBCalc_GetZreFromTau(this,CMB, tau) result(GetZreFromTau)
    use Reionization
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau
    type(CAMBParams)  P

    call this%CMBToCAMB(CMB, P)
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
        GetZreFromTau = Reion%GetZreFromTau(P, tau)
    class default
        call MpiStop('GetZreFromTau: Unknown reionization model')
    end select
    end function CAMBCalc_GetZreFromTau

    function CAMBCalc_CMBToTheta(this,CMB) result(CMBToTheta)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) CMBToTheta
    integer error

    call this%InitCAMB(CMB,error,.false.)
    CMBToTheta = this%CurrentState%CosmomcTheta()

    end function CAMBCalc_CMBToTheta


    real(mcp) function CAMBCalc_BAO_D_v(this, z)
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_BAO_D_v = this%CurrentState%BAO_D_v(z)

    end function CAMBCalc_BAO_D_v


    real(mcp) function CAMBCalc_AngularDiameterDistance(this, z)
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_AngularDiameterDistance = this%CurrentState%AngularDiameterDistance(z)

    end function CAMBCalc_AngularDiameterDistance

    real(mcp) function CAMBCalc_ComovingRadialDistance(this, z)
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_ComovingRadialDistance = this%CurrentState%ComovingRadialDistance(z)

    end function CAMBCalc_ComovingRadialDistance

    subroutine CAMBCalc_ComovingRadialDistanceArr(this, z, arr, n)
    class(CAMB_Calculator) :: this
    integer, intent(in) :: n
    real(mcp), intent(IN) :: z(n)
    real(mcp), intent(out) :: arr(n)
    !Note redshifts must be monotonically increasing

    call this%CurrentState%ComovingRadialDistanceArr(arr, z, n, 1d-4)

    end subroutine CAMBCalc_ComovingRadialDistanceArr

    real(mcp) function CAMBCalc_AngularDiameterDistance2(this, z1, z2)
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z1, z2

    CAMBCalc_AngularDiameterDistance2 = this%CurrentState%AngularDiameterDistance2(z1, z2)

    end function CAMBCalc_AngularDiameterDistance2

    real(mcp) function CAMBCalc_LuminosityDistance(this, z)
    !! distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_LuminosityDistance = this%CurrentState%LuminosityDistance(z)

    end function CAMBCalc_LuminosityDistance

    real(mcp) function CAMBCalc_Hofz(this, z)
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_Hofz = this%CurrentState%Hofz(z)

    end function CAMBCalc_Hofz


    subroutine CAMBCalc_InitCAMBParams(this,P)
    use lensing
    use DarkEnergyPPF
    use model
    use Reionization
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    integer zix
    !JD Changed P%Transfer%redshifts and P%Transfer%num_redshifts to
    !P%Transfer%PK_redshifts and P%Transfer%PK_num_redshifts respectively
    !for nonlinear lensing of CMB + LSS compatibility
    Threadnum =num_threads

    call CAMB_SetDefParams(P)
    if (this%dark_energy_model == 'ppf') then
        deallocate(P%DarkEnergy)
        allocate(TDarkEnergyPPF::P%DarkEnergy)
    else if (this%dark_energy_model /= 'fluid') then
        error stop 'Calculator_CAMB: unknown dark_energy_model'
    end if

    P%OutputNormalization = outNone

    !JD Modified to save computation time when only using MPK
    if(CosmoSettings%use_CMB) then
        P%WantScalars = .true.
        P%WantTensors = CosmoSettings%compute_tensors
        P%Max_l=CosmoSettings%lmax_computed_cl
        P%Max_eta_k=CosmoSettings%lmax_computed_cl*2
        P%Max_l_tensor=CosmoSettings%lmax_tensor
        P%Max_eta_k_tensor=CosmoSettings%lmax_tensor*5./2
    else
        P%WantCls = .false.
    end if

    P%WantTransfer = CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8
    P%Transfer%k_per_logint=0

    if (CosmoSettings%use_nonlinear) then
        P%NonLinear = NonLinear_pk
        P%Transfer%kmax = max(1.2_mcp,CosmoSettings%power_kmax)
    else
        P%Transfer%kmax = max(0.8_mcp,CosmoSettings%power_kmax)
    end if

    if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
        P%Transfer%high_precision=.true.
        P%Transfer%kmax=P%Transfer%kmax + 0.2
    end if
    P%Accuracy%AccuracyBoost = AccuracyLevel
    P%Accuracy%lAccuracyBoost = AccuracyLevel
    P%Accuracy%lSampleBoost = AccuracyLevel
    P%Accuracy%AccurateReionization = .true.

    P%Accuracy%AccurateBB = this%accurate_BB

    if (max_transfer_redshifts < CosmoSettings%num_power_redshifts) then
        stop 'Need to manually set max_transfer_redshifts larger in CAMB''s modules.f90'
    end if

    if (CosmoSettings%num_power_redshifts>1) then
        P%Transfer%PK_num_redshifts = CosmoSettings%num_power_redshifts
        do zix=1, CosmoSettings%num_power_redshifts
            !CAMB's ordering is from highest to lowest
            P%Transfer%PK_redshifts(zix) = CosmoSettings%power_redshifts(CosmoSettings%num_power_redshifts-zix+1)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts(1) = 0
    end if

    P%Num_Nu_Massive = 3
    P%Num_Nu_Massless = 0.046
    P%Accuracy%AccuratePolarization = CosmoSettings%num_cls/=1
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
        Reion%use_optical_depth = .false.
    class default
        call MpiStop('Reionization not TTanhReionization')
    end select

    if (CosmoSettings%CMB_Lensing) then
        P%DoLensing = .true.
        P%Max_l = CosmoSettings%lmax_computed_cl +100 + 50 !+50 in case accuracyBoost>1 and so odd l spacing
        P%Max_eta_k = P%Max_l*2
    end if

    P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl*AccuracyLevel,P%Max_eta_k)
    if (CosmoSettings%CMB_Lensing .and. (CosmoSettings%use_lensing_potential .or. CosmoSettings%use_nonlinear_lensing)) &
        P%Max_eta_k = max(P%Max_eta_k, 14000*AccuracyLevel)
    !k_etamax=18000 give c_phi_phi accurate to sub-percent at L=1000, <4% at L=2000
    !k_etamax=10000 is just < 1% at L<=500
    if (this%k_eta_max_scalar>0) then
        P%Max_eta_k = this%k_eta_max_scalar
    end if
    if (CosmoSettings%CMB_Lensing .and. CosmoSettings%use_nonlinear_lensing) then
        P%WantTransfer = .true.
        P%NonLinear = NonLinear_lens
        if(CosmoSettings%use_nonlinear) P%NonLinear = NonLinear_both
    end if
    lensing_includes_tensors = .false.

    P%Scalar_initial_condition = initial_vector
    allocate(P%InitialConditionVector(initial_nummodes),source=0._dl)
    P%InitialConditionVector(initial_adiabatic) = -1

    if (allocated(P%z_outputs)) deallocate(P%z_outputs)
    allocate(P%z_outputs, source = CosmoSettings%z_outputs)

    end subroutine CAMBCalc_InitCAMBParams

    !Mapping between array of power spectrum parameters and CAMB
    subroutine CAMBCalc_SetCAMBInitPower(this,P,CMB)
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    class(CMBParams) CMB

    select type (InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
        InitPower%pivot_scalar = CosmoSettings%pivot_k
        InitPower%pivot_tensor = CosmoSettings%tensor_pivot_k
        if (InitPower%pivot_tensor/=InitPower%pivot_scalar) InitPower%tensor_parameterization = tensor_param_rpivot
        InitPower%As = cl_norm*CMB%InitPower(As_index)
        InitPower%r = CMB%InitPower(amp_ratio_index)

        InitPower%ns = CMB%InitPower(ns_index)

        if (InitPower%r>0 .and. .not. CosmoSettings%compute_tensors) &
            call MpiStop('computing r>0 but compute_tensors=F')
        InitPower%nrun = CMB%InitPower(nrun_index)
        InitPower%nrunrun = CMB%InitPower(nrunrun_index)

        if (CosmoSettings%inflation_consistency) then
            if (CMB%InitPower(nt_index)/=0 .or. CMB%InitPower(ntrun_index)/=0) &
                & call MpiStop('Error: inflation_consistency but n_t not set to zero')
            ! first-order consistency relation
            !P%InitPower%At = - CMB%InitPower(amp_ratio_index)/8
            !next order consistency relation
            InitPower%At = - CMB%InitPower(amp_ratio_index)/8*(2-CMB%InitPower(ns_index) - CMB%InitPower(amp_ratio_index)/8)
            InitPower%ntrun = CMB%InitPower(amp_ratio_index)/8* &
                & (CMB%InitPower(amp_ratio_index)/8 + CMB%InitPower(ns_index) - 1)
            !note input n_T, nt run is ignored, so should be fixed
        else
            InitPower%At = CMB%InitPower(nt_index)
            InitPower%ntrun = CMB%InitPower(ntrun_index)
        end if
    class default
        stop 'CAMB_Calculator:Wrong initial power spectrum'
    end select

    end subroutine CAMBCalc_SetCAMBInitPower

    subroutine CAMBCalc_ReadParams(this,Ini)
    use NonLinear
    class(CAMB_Calculator) :: this
    class(TSettingIni) :: Ini


    call this%TCosmologyCalculator%ReadParams(Ini)
    this%calcName ='CAMB'

    this%CAMB_timing = Ini%Read_Logical('CAMB_timing',.false.)

    this%highL_theory_cl_template_file = Ini%ReadFilename('highL_theory_cl_template',DataDir,.true.)

    if (Ini%HasKey('highL_unlensed_cl_template')) then
        highL_unlensed_cl_template=  Ini%ReadFilename('highL_unlensed_cl_template')
    else
        highL_unlensed_cl_template = concat(LocalDir,'camb/fortran/',highL_unlensed_cl_template)
    end if

    this%k_eta_max_scalar = Ini%Read_Double('k_eta_max_scalar',-1._mcp)
    this%accurate_BB = Ini%Read_Logical('accurate_BB',.false.)

    this%halofit_version = Ini%Read_Int('halofit_version',halofit_default)

    this%dark_energy_model = Ini%Read_String_Default('dark_energy_model','ppf')

    end subroutine CAMBCalc_ReadParams


    subroutine CAMBCalc_InitForLikelihoods(this)
    !Called later after likelihoods etc loaded
    class(CAMB_Calculator) :: this

    if (CosmoSettings%use_CMB .and. CosmoSettings%lmax_computed_cl /= CosmoSettings%lmax) then
        if (CosmoSettings%compute_tensors .and. CosmoSettings%lmax_tensor > CosmoSettings%lmax_computed_cl) &
            & call MpiStop('lmax_tensor > lmax_computed_cl')
        call this%LoadFiducialHighLTemplate()
    end if

    call this%InitCAMBParams(this%CAMBP)
    select type(NL=>this%CAMBP%NonLinearModel)
    class is (THalofit)
        NL%halofit_version = this%halofit_version
    end select

    if (Feedback > 0 .and. MPIRank==0) then
        if(CosmoSettings%use_CMB) write(*,*) 'max_eta_k         = ', real(this%CAMBP%Max_eta_k)
        write(*,*) 'transfer kmax     = ', real(this%CAMBP%Transfer%kmax)
    end if

    this%CAMBP%WantTensors = CosmoSettings%compute_tensors

    end subroutine CAMBCalc_InitForLikelihoods


    subroutine CAMBCalc_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    use Recombination
    use Reionization
    class(CAMB_Calculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.
    call ReadValues%Add('Compiled_CAMB_version', version)
    call ReadValues%Add('Compiled_Equations', Eqns_name)

    end subroutine CAMBCalc_VersionTraceOutput

    subroutine CAMBcalc_SetCurrentPoint(this, Params)
    class(CAMB_Calculator), target :: this
    class(TCalculationAtParamPoint), target :: Params

    if (.not. associated(Params%Info)) then
        allocate(Params%Info, source=this%DefaultInstance)
    end if

    select type(Info=>Params%Info)
    class is (CAMBTransferCache)
        this%CurrentState => Info%State
    end select

    end subroutine CAMBcalc_SetCurrentPoint

    subroutine LoadFiducialHighLTemplate(this)
    class(CAMB_Calculator) :: this
    !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
    integer L,  status
    real(mcp) array(4)
    Type(TTextFile) :: F

    allocate(this%highL_lensedCL_template(2:CosmoSettings%lmax, 4))
    call F%Open(this%highL_theory_cl_template_file)
    do
        read(F%unit,*, iostat=status) L , array
        if (status/=0 .or. L>CosmoSettings%lmax) exit
        if (L>=2) this%highL_lensedCL_template(L,:) = array
    end do
    call F%Close()

    if (this%highL_lensedCL_template(2,1) < 100) &
        call MpiStop('highL_theory_cl_template must be in muK^2')

    if (L<CosmoSettings%lmax) call MpiStop('highL_theory_cl_template does not go to lmax')

    end subroutine LoadFiducialHighLTemplate

    end module
