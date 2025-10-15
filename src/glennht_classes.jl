module GlennHTClasses

# ========== Enums ==========
@enum BoundaryConditionType begin
    InletBCType = 0
    OutletBCType = 1
    SymmetryOrSlipBCType = 2
    WallBCType = 3
    GIFBCType = 4
end

@enum InletBC_Subtype  Normal=0  AngleSpecified=1  AngleAndProfileSpecified=2
@enum InletBC_Direction Something=1  Annular=2  Cascade=3
@enum OutletSubtype  UniformPressure=0  AveragePressure=1  UnsteadyFlow=2
@enum SymmetricSlipSubtype  Symmetry=0  Slip=1
@enum WallSubtype  SpecifiedWallHeatFlux=0  SpecifiedWallTemperature=1  BCWall_Subtype_Conjugate=3
@enum GIFCoordinate  Cartesian=0  Polar=1
@enum GIFType  Conjugate=1  StraightInterpolation=2
@enum GIFOrder  ZeroOrder=0  LinearOrder=1  CubicOrder=3
@enum TbModelType  NO_TURB_MODEL=0  K_OMEGA_TURB_MODEL=1  K_EPS_TURB_MODEL=2  ARSM_TURB_MODEL=3  RSM_TURB_MODEL=4  SST_TURB_MODEL=5  LES_TURB_MODEL=7  WL_KOMEGA_TURB_MODEL=11  K_OMEGA_GAMMA_TURB_MODEL=12

# ========== Boundary conditions ==========
abstract type AbstractBoundaryCondition end

mutable struct InletBC <: AbstractBoundaryCondition
    BCType::BoundaryConditionType
    SurfaceID::Int
    Name::String
    IsPostProcessing::Bool
    IsCalculateMassFlow::Bool
    ToggleProcessSurface::Bool

    inlet_subType::InletBC_Subtype
    inlet_ref_Mach_Nr::Float64

    T0_const::Union{Nothing,Float64}
    P0_const::Union{Nothing,Float64}
    P0_const_unit::String
    Tu_const::Union{Nothing,Float64}
    Ts_const::Union{Nothing,Float64}
    ang1_const::Union{Nothing,Float64}
    bet1_const::Union{Nothing,Float64}

    annular_inlet::Bool
    deltah::Union{Nothing,Float64}
    deltat::Union{Nothing,Float64}
    twall_hub::Union{Nothing,Float64}
    twall_case::Union{Nothing,Float64}

    have_inlet_prof::Bool
    filen_inlet_profile::Union{Nothing,String}
    direction::Union{Nothing,InletBC_Direction}

    function InletBC(; SurfaceID::Int, Name::String="",
                     IsPostProcessing::Bool=false, IsCalculateMassFlow::Bool=false, ToggleProcessSurface::Bool=false,
                     inlet_subType::InletBC_Subtype=Normal, inlet_ref_Mach_Nr::Float64=0.2,
                     T0_const=nothing, P0_const=nothing, P0_const_unit::String="Pa", Tu_const=nothing, Ts_const=nothing,
                     ang1_const=nothing, bet1_const=nothing,
                     annular_inlet::Bool=false, deltah=nothing, deltat=nothing, twall_hub=nothing, twall_case=nothing,
                     have_inlet_prof::Bool=false, filen_inlet_profile=nothing, direction=nothing)
        new(InletBCType, SurfaceID, Name, IsPostProcessing, IsCalculateMassFlow, ToggleProcessSurface,
            inlet_subType, inlet_ref_Mach_Nr,
            T0_const, P0_const, P0_const_unit, Tu_const, Ts_const, ang1_const, bet1_const,
            annular_inlet, deltah, deltat, twall_hub, twall_case,
            have_inlet_prof, filen_inlet_profile, direction)
    end
end

mutable struct OutletBC <: AbstractBoundaryCondition
    BCType::BoundaryConditionType
    SurfaceID::Int
    Name::String
    IsPostProcessing::Bool
    IsCalculateMassFlow::Bool
    ToggleProcessSurface::Bool

    outlet_subType::Union{Nothing,OutletSubtype}
    extrapolation_order::Union{Nothing,Int}
    Pback_extrapolate_profile::Bool
    Pback_const::Union{Nothing,Float64}
    Pback_const_unit::String
    have_Pback_prof::Bool
    annular_outlet::Bool
    approx_Mach_out::Union{Nothing,Float64}
    filen_Pback_prof::Union{Nothing,String}

    mult_for_full_ring::Union{Nothing,Int}
    p_over_pt_ratio_limit::Union{Nothing,Float64}

    function OutletBC(; SurfaceID::Int, Name::String="",
                      IsPostProcessing::Bool=false, IsCalculateMassFlow::Bool=false, ToggleProcessSurface::Bool=false,
                      outlet_subType=nothing, extrapolation_order=nothing, Pback_extrapolate_profile::Bool=false,
                      Pback_const=nothing, Pback_const_unit::String="Pa", have_Pback_prof::Bool=false,
                      annular_outlet::Bool=false, approx_Mach_out=nothing, filen_Pback_prof=nothing,
                      mult_for_full_ring=nothing, p_over_pt_ratio_limit=nothing)
        new(OutletBCType, SurfaceID, Name, IsPostProcessing, IsCalculateMassFlow, ToggleProcessSurface,
            outlet_subType, extrapolation_order, Pback_extrapolate_profile,
            Pback_const, Pback_const_unit, have_Pback_prof, annular_outlet, approx_Mach_out, filen_Pback_prof,
            mult_for_full_ring, p_over_pt_ratio_limit)
    end
end

mutable struct SymmetricSlipBC <: AbstractBoundaryCondition
    BCType::BoundaryConditionType
    SurfaceID::Int
    Name::String
    IsPostProcessing::Bool
    IsCalculateMassFlow::Bool
    ToggleProcessSurface::Bool

    slip_subType::Union{Nothing,SymmetricSlipSubtype}
    slip_omega::Union{Nothing,Float64}

    function SymmetricSlipBC(; SurfaceID::Int, Name::String="",
                             IsPostProcessing::Bool=false, IsCalculateMassFlow::Bool=false, ToggleProcessSurface::Bool=false,
                             slip_subType=nothing, slip_omega=nothing)
        new(SymmetryOrSlipBCType, SurfaceID, Name, IsPostProcessing, IsCalculateMassFlow, ToggleProcessSurface,
            slip_subType, slip_omega)
    end
end

mutable struct WallBC <: AbstractBoundaryCondition
    BCType::BoundaryConditionType
    SurfaceID::Int
    Name::String
    IsPostProcessing::Bool
    IsCalculateMassFlow::Bool
    ToggleProcessSurface::Bool

    wall_subType::Union{Nothing,Int}
    Twall_const::Union{Nothing,Float64}
    have_Twall_prof::Bool
    filen_Twall_prof::Union{Nothing,String}
    Qwall_const::Union{Nothing,Float64}
    have_Qwall_prof::Bool
    filen_Qwall_prof::Union{Nothing,String}
    BEM_coupled_surf::Bool
    Nr_wall_segments::Union{Nothing,Int}
    segment_Omega::Union{Nothing,Float64}
    segment_xMin::Union{Nothing,Float64}

    function WallBC(; SurfaceID::Int, Name::String="",
                    IsPostProcessing::Bool=false, IsCalculateMassFlow::Bool=false, ToggleProcessSurface::Bool=false,
                    wall_subType=nothing, Twall_const=nothing, have_Twall_prof::Bool=false, filen_Twall_prof=nothing,
                    Qwall_const::Union{Nothing,Float64}=0.0, have_Qwall_prof::Bool=false, filen_Qwall_prof=nothing,
                    BEM_coupled_surf::Bool=false, Nr_wall_segments=nothing, segment_Omega=nothing, segment_xMin=nothing)
        new(WallBCType, SurfaceID, Name, IsPostProcessing, IsCalculateMassFlow, ToggleProcessSurface,
            wall_subType, Twall_const, have_Twall_prof, filen_Twall_prof, Qwall_const, have_Qwall_prof, filen_Qwall_prof,
            BEM_coupled_surf, Nr_wall_segments, segment_Omega, segment_xMin)
    end
end

mutable struct GIF <: AbstractBoundaryCondition
    BCType::BoundaryConditionType
    SurfaceID::Int             # not used directly; pair below is authoritative
    Name::String
    IsPostProcessing::Bool
    IsCalculateMassFlow::Bool
    ToggleProcessSurface::Bool

    GIFSurface1::Int
    GIFSurface2::Int
    Name1::String
    Name2::String
    Coordinates::GIFCoordinate
    Type::GIFType
    Order::GIFOrder

    function GIF(; GIFSurface1::Int, GIFSurface2::Int, Name1::String="", Name2::String="",
                  Coordinates::GIFCoordinate=Cartesian, Type::GIFType=Conjugate, Order::GIFOrder=LinearOrder)
        new(GIFBCType, 0, "", false, false, false,
            GIFSurface1, GIFSurface2, Name1, Name2, Coordinates, Type, Order)
    end
end

mutable struct BCGroup
    Inlets::Vector{InletBC}
    Outlets::Vector{OutletBC}
    SymmetricSlips::Vector{SymmetricSlipBC}
    Walls::Vector{WallBC}
    function BCGroup(; Inlets=InletBC[], Outlets=OutletBC[], SymmetricSlips=SymmetricSlipBC[], Walls=WallBC[])
        new(Inlets, Outlets, SymmetricSlips, Walls)
    end
end

# ========== Job-related structures (minimal, mirroring Python names) ==========
mutable struct JobFiles
    DcmpFILE::Union{Nothing,String}
    ConnFILE::Union{Nothing,String}
    BCSpecFILE::Union{Nothing,String}
    GridFile::Union{Nothing,String}
    GridFileFormat::Union{Nothing,String}
    Plot3DFileFormat::Union{Nothing,String}
    SolnInFile::Union{Nothing,String}
    SolnInFileFormat::Union{Nothing,String}
    SolnOutFile::Union{Nothing,String}
    SolnOutFileFormat::Union{Nothing,String}
    residFILE::Union{Nothing,String}
    residFILE2::Union{Nothing,String}
    function JobFiles()
        new("ddcmp.dat", "connectivity.ght_conn", "boundary_conditions.bcs", "mesh.xyz", "formatted",
            "formatted", "In.soln", "unformatted", "Out.soln", "unformatted", "his.subs", "his.nosubs")
    end
end

mutable struct JobControl
    mRunLevel::Int
    LUNout::Int
    RestartSoln::Bool
    SaveSoln::Bool
    SaveTransientSoln::Bool
    VerboseScreenOutput::Bool
    JobControl() = new(0, 6, false, true, false, true)
end

mutable struct TurbModelInput
    TbModelType_::TbModelType
    PRNS_ResolutionParameter::Float64
    TurbModelInput() = new(K_OMEGA_TURB_MODEL, 1.0)
end

mutable struct Plot3DParameters
    Plot3DParameterSet::String
    Plot3DParameters() = new("Standard")
end

mutable struct InitialCond
    P0::Float64; T0::Float64; Minit::Float64; alfa::Float64; beta::Float64
    Tu::Float64; Ts::Float64; T0_solid::Float64; annular_init::Float64
    InitialCond() = new(1.0,1.0,0.0,0.0,0.0,0.05,0.05,1.0,0.0)
end

mutable struct TimeStpControl
    UnsteadyFlow::Bool; FullyImplicitDiscr::Bool; EulerBackward::Bool; CranckNicholson::Bool
    BlendedTime::Bool; TransientPlot3dFiles::Bool; have_previous_step::Bool; UseLowMPrecond::Bool
    pcMinRefMach::Float64; dissRefMach::Float64; Implicitness::Float64; dt_unst::Float64
    CFLn::Float64; CFLr::Float64; cst::Float64; cst_solid::Float64
    convergTolerance::Float64; nTimeSteps::Int; maxPseudoSteps::Int; maxPseudoSteps_solid::Int
    fully_coupled_solid::Bool; ReinitializeTime::Bool; ResetTimeTo::Float64
    nTransBegin::Int; nfiles::Int; ninterval::Int; nHL::Int
    time_avg_start::Union{Nothing,Float64}; time_avg_end::Union{Nothing,Float64}; restart_timeAvg::Bool
    TimeStpControl() = new(false,false,false,false,true,false,false,false,
                           1e-3,1e-3,1.0,1e-8, 0.25,0.125,3.5,7.0,
                           5e-5,50,50,100, false,false,0.0, 0,10,10,1, nothing,nothing,false)
end

mutable struct SPDSchemeControl
    NS_Central::String; TB2_Upwind1::String; NS_Upwind2::String
    ScalrCoeff_ArtDiss::Bool; useSecDiffArtDiss::Bool; useFrthDiffArtDiss::Bool
    rk2::String; rk4::String
    NS_Upwind1::Bool; use_AUSM_Chima::Bool; use_AUSM_Liou_hTot::Bool; TB2_Central::Bool
    TBRSM_Central::Bool; TBRSM_Upwind1::Bool; constArtDiss::Bool; scalarArtDiss::Bool
    MatrxCoeff_ArtDiss::Bool; secDiffArtDiss::Bool; matrixArtDiss::Bool; frthDiffArtDiss::Bool
    MachCutOff::Float64; ivanAlbada::Int
    SPDSchemeControl() = new("4*T","4*T","4*F", true,true,true, "4*0.12500","4*0.032",
                             false,false,false,false, false,true,false,true, false,true,false,true, 0.1,1)
end

mutable struct RKSchemeControl
    nStages::Int; RKCoeff::String; compute_pdiff_in_stage::String
    compute_adiss_in_stage::String; export_import_after_stage::String
    use_implicit_residual_smoothing::String; irs_neqs::Union{Nothing,Int}
    irs_use_GS::Bool; n_GS_iterations::Union{Nothing,Int}; n_GlobalSweeps::Union{Nothing,Int}
    RKSchemeControl() = new(4, "0.25,0.3333333,0.5,1.,6*0", "T,T,T,T,6*F",
                            "T,T,T,T,6*F", "T,T,T,T,6*F", ".T.", 1, true, 3, 1)
end

mutable struct MGSchemeControl
    FinestLevel::Int; CoarsestLevel::Int; pre_mg_sweeps::Int; mg_sweeps::Int; post_mg_sweeps::Int
    SVFinestLevel::Int; SVCoarsestLevel::Int; SVpre_mg_sweeps::Int; SVmg_sweeps::Int; SVpost_mg_sweeps::Int
    MGSchemeControl() = new(0,0,1,0,0, 0,0,1,0,0)
end

mutable struct GasPropertiesInput
    UseDryAir::Bool; Use_const_Cp::Bool; Use_const_trProp::Bool; Use_RefT::Bool
    RefT_Properties::Float64; const_cp::Float64; const_visc::Float64; const_kth::Float64
    Use_specialGas::Bool; SpecialGasMW::Float64
    GasPropertiesInput() = new(true,true,false,true, -1.0e99,-1.0e99,-1.0e99,-1.0e99, false,-1.0e99)
end

mutable struct ReferenceCond
    useDimensionalVariables::Bool; refLen::Float64; refP0::Float64; refT0::Float64
    refRho0::Float64; refVel::Float64; refVisc::Float64; refCond::Float64; refCp::Float64
    MolW::Float64; RgasUnv::Float64; Rgas::Float64; gamma::Float64; Re::Float64; Pr::Float64
    ndVisc::Float64; ndCond::Float64; Omegab::Float64; ReScalingFactor::Float64
    rho_solid::Union{Nothing,Float64}; cond_solid::Union{Nothing,Float64}; Csp_solid::Union{Nothing,Float64}
end

mutable struct ReferenceCondFull
    reflen::Union{Nothing,Float64}; refP0::Union{Nothing,Float64}; refT0::Union{Nothing,Float64}
    refrho0::Union{Nothing,Float64}; refVel::Union{Nothing,Float64}; refvisc::Union{Nothing,Float64}
    refcond::Union{Nothing,Float64}; refCp::Union{Nothing,Float64}; MolW::Union{Nothing,Float64}
    RgasUnv::Union{Nothing,Float64}; Rgas::Union{Nothing,Float64}; gamma::Union{Nothing,Float64}
    Re::Union{Nothing,Float64}; Pr::Union{Nothing,Float64}; ndVisc::Union{Nothing,Float64}
    ndCond::Union{Nothing,Float64}; Omegab::Union{Nothing,Float64}; ReScalingFactor::Union{Nothing,Float64}
    rho_solid::Union{Nothing,Float64}; cond_solid::Union{Nothing,Float64}; csp_solid::Union{Nothing,Float64}
    ReferenceCondFull() = new(nothing,nothing,nothing, nothing,nothing,nothing, nothing,nothing,nothing,
                              nothing,nothing,nothing, nothing,nothing,nothing, nothing,nothing,nothing,
                              nothing,nothing,nothing)
end

mutable struct Job
    JobFiles::JobFiles
    JobControl::JobControl
    TurbModelInput::TurbModelInput
    Plot3DParameters::Plot3DParameters
    InitialCond::InitialCond
    TimeStpControl::TimeStpControl
    SPDSchemeControl::SPDSchemeControl
    RKSchemeControl::RKSchemeControl
    MGSchemeControl::MGSchemeControl
    GasPropertiesInput::GasPropertiesInput
    ReferenceCondFull::ReferenceCondFull
    ReferenceCond::Union{Nothing,ReferenceCond}
    Job() = new(JobFiles(), JobControl(), TurbModelInput(), Plot3DParameters(), InitialCond(),
                TimeStpControl(), SPDSchemeControl(), RKSchemeControl(), MGSchemeControl(),
                GasPropertiesInput(), ReferenceCondFull(), nothing)
end

end # module
