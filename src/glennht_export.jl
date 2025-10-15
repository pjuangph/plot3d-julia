import JSON                      # if you serialize a sidecar; or drop if not used
import Printf: @printf

import .GlennHTClasses:
    Job, BCGroup, GIF,
    AbstractBoundaryCondition,
    ReferenceCond, ReferenceCondFull

export to_pa, ideal_R, mach_from_p0_over_p, T_from_T0, mu_suth, a_sound, cp_from_gamma_R,
       pr_and_k, rpm_to_omegab, populate_reference_from_inputs,
       export_to_boundary_condition, export_to_job_file, summarize_contiguous,
       export_to_glennht_conn, ensure_extension

const _AIR_MOLW_DEFAULT = 28.964
const _R_UNIV = 8314.4126
const _SUTH_mu0 = 1.716e-5
const _SUTH_T0  = 273.15
const _SUTH_S   = 110.4

const _PRESSURE_TO_PA = Dict(
    "pa"=>1.0,"kpa"=>1e3,"mpa"=>1e6,"bar"=>1e5,"mbar"=>1e2,
    "atm"=>101325.0,"psi"=>6894.757293168,
)

to_pa(value::Union{Nothing,Real}, unit::Union{Nothing,AbstractString}) =
    value === nothing ? nothing : (let u = lowercase(String(unit === nothing ? "Pa" : unit));
        haskey(_PRESSURE_TO_PA, u) || error("Unknown pressure unit '$unit'")
        float(value) * _PRESSURE_TO_PA[u]
    end)

ideal_R(molw::Real) = _R_UNIV/float(molw)
mach_from_p0_over_p(p0_over_p::Union{Nothing,Real}, gamma::Real) =
    (p0_over_p === nothing || p0_over_p <= 1.0) ? 0.0 :
    sqrt(max(0.0, 2.0 * (p0_over_p^((gamma-1)/gamma) - 1.0) / (gamma-1.0)))
T_from_T0(T0::Real, M::Real, γ::Real) = T0/(1.0 + 0.5*(γ-1.0)*M*M)
mu_suth(T::Real) = _SUTH_mu0 * (T/_SUTH_T0)^1.5 * (_SUTH_T0 + _SUTH_S)/(T + _SUTH_S)
a_sound(T::Real, γ::Real, R::Real) = sqrt(γ*R*T)
cp_from_gamma_R(γ::Real, R::Real) = γ*R/(γ-1.0)
function pr_and_k(mu::Real, cp::Real, pr_given::Union{Nothing,Real}, k_given::Union{Nothing,Real})
    if pr_given !== nothing && k_given === nothing
        return (float(pr_given), (mu*cp)/float(pr_given))
    end
    if pr_given === nothing && k_given !== nothing && k_given != 0
        return ((mu*cp)/float(k_given), float(k_given))
    end
    pr_default = 0.706
    return (pr_default, (mu*cp)/pr_default)
end
rpm_to_omegab(rpm::Union{Nothing,Real}) = (rpm === nothing || rpm == 0) ? 0.0 : 2π * (float(rpm)/60.0)

ensure_extension(path::AbstractString, ext::AbstractString=".txt") =
    endswith(lowercase(path), lowercase(ext)) ? path : string(path[begin:end], ext)

# ---------------- Reference builder ----------------
function populate_reference_from_inputs(job::GlennHTClasses.Job;
    inlet_total_P0_Pa::Union{Nothing,Real},
    reference_static_P_Pa::Union{Nothing,Real},
    inlet_T0_K::Union{Nothing,Real},
    refLen_m::Union{Nothing,Real},
    MolW::Union{Nothing,Real}=nothing,
    gamma::Union{Nothing,Real}=nothing,
    Pr::Union{Nothing,Real}=nothing,
    k_override_WmK::Union{Nothing,Real}=nothing,
    rpm::Union{Nothing,Real}=nothing,
    rho_solid::Union{Nothing,Real}=nothing,
)
    rcfull = job.ReferenceCondFull
    MolW  = MolW  === nothing || MolW  == 0 ? _AIR_MOLW_DEFAULT : float(MolW)
    gamma = gamma === nothing || gamma == 0 ? 1.4 : float(gamma)
    R     = ideal_R(MolW)
    T0    = inlet_T0_K === nothing || inlet_T0_K == 0 ? (rcfull.refT0 === nothing ? 300.0 : rcfull.refT0) : float(inlet_T0_K)

    rcfull.refP0 = inlet_total_P0_Pa === nothing || inlet_total_P0_Pa == 0 ? (rcfull.refP0 === nothing ? 101325.0 : rcfull.refP0) : float(inlet_total_P0_Pa)
    rcfull.reflen = refLen_m === nothing || refLen_m == 0 ? (rcfull.reflen === nothing ? 1.0 : rcfull.reflen) : float(refLen_m)

    M = (reference_static_P_Pa !== nothing && rcfull.refP0 !== nothing) ? mach_from_p0_over_p(rcfull.refP0 / float(reference_static_P_Pa), gamma) : 0.0
    T_static = T_from_T0(T0, M, gamma)
    a = a_sound(T_static, gamma, R)
    refVel = M * a
    mu = mu_suth(T_static)
    Pref = reference_static_P_Pa === nothing ? rcfull.refP0 : float(reference_static_P_Pa)
    rho = Pref / (R*T_static)

    cp = cp_from_gamma_R(gamma, R)
    (Pr_val, k_val) = pr_and_k(mu, cp, Pr === nothing ? nothing : float(Pr),
                                k_override_WmK === nothing ? nothing : float(k_override_WmK))

    Re = (refVel != 0 && mu != 0 && rcfull.reflen !== nothing) ? (rho * refVel * float(rcfull.reflen) / mu) : 0.0

    rcfull.refT0   = T0
    rcfull.refrho0 = rho
    rcfull.refVel  = refVel
    rcfull.refvisc = mu
    rcfull.refcond = k_val
    rcfull.refCp   = cp
    rcfull.MolW    = MolW
    rcfull.RgasUnv = _R_UNIV
    rcfull.Rgas    = R
    rcfull.gamma   = gamma
    rcfull.Re      = Re
    rcfull.Pr      = Pr_val
    rcfull.ndVisc  = 1.0
    rcfull.ndCond  = 1.0
    rcfull.Omegab  = rpm_to_omegab(rpm)
    rcfull.ReScalingFactor = 1.0
    rcfull.rho_solid = rho_solid === nothing ? nothing : float(rho_solid)
    rcfull.cond_solid = 20.0
    rcfull.csp_solid  = 896.0

    job.ReferenceCond = GlennHTClasses.ReferenceCond(false,
        rcfull.reflen, rcfull.refP0, rcfull.refT0, rcfull.refrho0, rcfull.refVel, rcfull.refvisc,
        rcfull.refcond, rcfull.refCp, rcfull.MolW, rcfull.RgasUnv, rcfull.Rgas, rcfull.gamma,
        rcfull.Re, rcfull.Pr, 1.0, 1.0, rcfull.Omegab, 1.0,
        rcfull.rho_solid, rcfull.cond_solid, rcfull.csp_solid)
    return nothing
end

# ---------------- Formatting helpers ----------------
_fmt_bool(v::Bool) = v ? ".TRUE." : ".FALSE."
function _fmt_value(v)
    if v === nothing; return nothing; end
    if v isa Bool; return _fmt_bool(v); end
    if v isa Real; return repr(float(v)); end
    if v isa AbstractString; return "'$v'"; end
    if v isa Tuple || v isa Vector
        return join(filter(!isnothing, (_fmt_value(x) for x in v)), ",")
    end
    return "'$(string(v))'"
end

function _iter_fields(obj)
    ks = fieldnames(typeof(obj))
    ((String(k), getfield(obj, k)) for k in ks)
end

function _export_namelist_block(header::AbstractString, obj; exclude_names=Set{String}())
    pairs = String[]
    for (k,v) in _iter_fields(obj)
        (v === nothing || (k in exclude_names) || endswith(k, "_unit")) && continue
        fv = _fmt_value(v); fv === nothing && continue
        push!(pairs, "$k=$fv")
    end
    inner = join(pairs, ", ")
    return " &$header\n$inner\n &END\n"
end

function _write_bsurf_spec(io::IO, bc::AbstractBoundaryCondition)
    line = " &BSurf_Spec\n" *
           "BSurfID=$(bc.SurfaceID), BCType=$(Int(bc.BCType)), BSurfName='$(bc.Name)'"
    if bc.IsPostProcessing; line *= ", BRefCond=T"; end
    if bc.IsCalculateMassFlow || bc.ToggleProcessSurface; line *= ", BCalc=T"; end
    write(io, line * "\n &END\n\n")
end

# ---------------- GIF + Volume Zones ----------------
function _write_gif_from_dict(io::IO, gdict::Dict{String,Any})
    sid1 = Int(get(gdict, "a", 0))
    sid2 = Int(get(gdict, "b", 0))
    name1 = "surface $sid1"; name2 = "surface $sid2"; bctype = 4
    write(io, " &BSurf_Spec\nBSurfID=$sid1, BCType=$bctype, BSurfName='$name1'\n &END\n\n")
    write(io, " &BSurf_Spec\nBSurfID=$sid2, BCType=$bctype, BSurfName='$name2'\n &END\n\n")
    write(io, " &GIF_Spec\nSurfID_1=$sid1, SurfID2=$sid2\n &END\n\n")
end

function _write_vzconditions(io::IO, vz::Dict{String,Any})
    vzid = Int(get(vz, "contiguous_id", 0))
    ztype = lowercase(String(get(vz, "zone_type", "fluid")))
    vztype = ztype == "fluid" ? 1 : 2
    if vztype == 1
        write(io,
            " &VZConditions\n" *
            "VZid=$vzid, VZtype=1, OmegaVZ=0., VZMaterialName=Air,\n" *
            "Fluid_Tref_prop=0., Fluid_k_Tref=285., Fluid_amu_Tref=285., Fluid_expnt=.7,UseDryAir=.TRUE.,\n" *
            "!Fluid_cp=1002., Fluid_Pr=.7, Fluid_MW=28.964\n" *
            " &END\n\n"
        )
    else
        write(io,
            " &VZConditions\n" *
            "VZid=$vzid, VZtype=2, OmegaVZ=0., VZMaterialName=CMC,\n" *
            "Solid_Tref_prop=285., Solid_rho_Tref=2707. , Solid_condN_Tref=6.5, Solid_condT_Tref=6.5, Solid_condA_Tref=6.5,\n" *
            "Solid_Csp_Tref=896.\n" *
            " &END\n\n"
        )
    end
end

# ---------------- BCS export ----------------
function export_to_boundary_condition(file_path_to_write::AbstractString,
                                      job_settings::GlennHTClasses.Job,
                                      bc_group::GlennHTClasses.BCGroup,
                                      gif_pairs::Vector{Union{GIF,Dict{String,Any}}},     # allow dicts or GIF objects
                                      volume_zones::Vector{Dict{String,Any}})
    path = ensure_extension(file_path_to_write, ".bcs")
    open(path, "w") do io
        # Robust reference defaults using first inlet if missing
        ref = job_settings.ReferenceCondFull
        first_inlet = isempty(bc_group.Inlets) ? nothing : bc_group.Inlets[1]
        if ref.reflen === nothing || ref.reflen == 0; ref.reflen = 1.0; end
        if (ref.refP0 === nothing || ref.refP0 == 0) && first_inlet !== nothing && first_inlet.P0_const !== nothing
            phys = to_pa(first_inlet.P0_const, first_inlet.P0_const_unit)
            ref.refP0 = phys
        end
        if (ref.refT0 === nothing || ref.refT0 == 0) && first_inlet !== nothing && first_inlet.T0_const !== nothing
            ref.refT0 = first_inlet.T0_const
        end

        # Normalize inlet values by references
        for inlet in bc_group.Inlets
            if inlet.P0_const !== nothing && ref.refP0 !== nothing
                phys = to_pa(inlet.P0_const, inlet.P0_const_unit)
                inlet.P0_const = phys/ref.refP0
            end
            if inlet.T0_const !== nothing && ref.refT0 !== nothing
                inlet.T0_const = inlet.T0_const/ref.refT0
            end
            if inlet.twall_hub !== nothing && ref.refT0 !== nothing
                inlet.twall_hub = inlet.twall_hub/ref.refT0
            end
            if inlet.twall_case !== nothing && ref.refT0 !== nothing
                inlet.twall_case = inlet.twall_case/ref.refT0
            end
            if inlet.Ts_const !== nothing && ref.reflen !== nothing
                inlet.Ts_const = inlet.Ts_const/ref.reflen
            end
            _write_bsurf_spec(io, inlet)
        end

        # Outlets
        for outlet in bc_group.Outlets
            if outlet.Pback_const !== nothing && ref.refP0 !== nothing
                phys = to_pa(outlet.Pback_const, outlet.Pback_const_unit)
                outlet.Pback_const = phys/ref.refP0
            end
            _write_bsurf_spec(io, outlet)
        end

        # Slips / Walls
        for s in bc_group.SymmetricSlips; _write_bsurf_spec(io, s); end
        for w in bc_group.Walls; _write_bsurf_spec(io, w); end

        # GIFs
        for pair in gif_pairs
            if pair isa Dict
                _write_gif_from_dict(io, pair)
            else
                p = pair::GlennHTClasses.GIF
                write(io, " &BSurf_Spec\nBSurfID=$(p.GIFSurface1), BCType=$(Int(p.BCType)), BSurfName='$(p.Name1)'\n &END\n\n")
                write(io, " &BSurf_Spec\nBSurfID=$(p.GIFSurface2), BCType=$(Int(p.BCType)), BSurfName='$(p.Name2)'\n &END\n\n")
                write(io, " &GIF_Spec\nSurfID_1=$(p.GIFSurface1), SurfID2=$(p.GIFSurface2)\n &END\n\n")
            end
        end

        # Unique volume zones by contiguous_id
        seen = Set{Int}()
        for vz in volume_zones
            cid = Int(vz["contiguous_id"])
            cid in seen && continue
            push!(seen, cid)
            _write_vzconditions(io, vz)
        end

        # Detailed BC blocks (only unique subtype representatives)
        function first_by_subtype(objs, getsub)
            seen = Set{Any}(); out = Any[]
            for o in objs
                st = getsub(o)
                (st === nothing) && continue
                st in seen && continue
                push!(seen, st); push!(out, o)
            end
            out
        end
        excl = Set(["Name","SurfaceID","BCType"])
        for inlet in first_by_subtype(bc_group.Inlets, o->o.inlet_subType)
            write(io, _export_namelist_block("INLET_BC", inlet; exclude_names=excl)); write(io,"\n")
        end
        for outlet in first_by_subtype(bc_group.Outlets, o->o.outlet_subType)
            write(io, _export_namelist_block("OUTLET_BC", outlet; exclude_names=excl)); write(io,"\n")
        end
        for slip in first_by_subtype(bc_group.SymmetricSlips, o->o.slip_subType)
            write(io, _export_namelist_block("SLIP_BC", slip; exclude_names=excl)); write(io,"\n")
        end
        for wall in first_by_subtype(bc_group.Walls, o->o.wall_subType)
            write(io, _export_namelist_block("WALL_BC", wall; exclude_names=excl)); write(io,"\n")
        end
    end
    return nothing
end

# ---------------- Job export ----------------
function export_to_job_file(job::GlennHTClasses.Job, file_path_to_write::AbstractString;
                            title::Union{Nothing,AbstractString}=nothing,
                            exec_serial::AbstractString="GlennHT.serial",
                            exec_mpi::Union{Nothing,AbstractString}=nothing)
    open(file_path_to_write, "w") do io
        if title !== nothing
            write(io, " &Title\n TheTitle=\"$title\"\n &end\n\n")
        end
        if job.ReferenceCond === nothing && job.ReferenceCondFull !== nothing
            try
                populate_reference_from_inputs(job;
                    inlet_total_P0_Pa = job.ReferenceCondFull.refP0,
                    reference_static_P_Pa = nothing,
                    inlet_T0_K = job.ReferenceCondFull.refT0,
                    refLen_m   = job.ReferenceCondFull.reflen)
            catch
            end
        end
        write(io, _export_namelist_block("JobFiles", job.JobFiles)); write(io,"\n")
        write(io, _export_namelist_block("JobControl", job.JobControl)); write(io,"\n")
        write(io, _export_namelist_block("TurbModelInput", job.TurbModelInput)); write(io,"\n")
        write(io, _export_namelist_block("Plot3DParameters", job.Plot3DParameters)); write(io,"\n")
        write(io, _export_namelist_block("InitialCond", job.InitialCond)); write(io,"\n")
        write(io, _export_namelist_block("TimeStpControl", job.TimeStpControl)); write(io,"\n")
        write(io, _export_namelist_block("SPDSchemeControl", job.SPDSchemeControl)); write(io,"\n")
        write(io, _export_namelist_block("RKSchemeControl", job.RKSchemeControl)); write(io,"\n")
        write(io, _export_namelist_block("MGSchemeControl", job.MGSchemeControl)); write(io,"\n")
        write(io, _export_namelist_block("GasPropertiesInput", job.GasPropertiesInput)); write(io,"\n")
        if job.ReferenceCond !== nothing
            write(io, _export_namelist_block("ReferenceCond", job.ReferenceCond)); write(io,"\n")
        end
        if job.ReferenceCondFull !== nothing
            write(io, _export_namelist_block("ReferenceCondFull", job.ReferenceCondFull)); write(io,"\n")
        end
        write(io, "execFILE=\"$exec_serial\"\n")
        if exec_mpi !== nothing
            write(io, "execFILE=\"$exec_mpi\"\n")
        end
    end
    return nothing
end

# ---------------- Summaries + GlennHT connectivity export ----------------
function summarize_contiguous(records::Vector{Dict{String,Any}})
    id_to_types = Dict{Int, Set{String}}()
    for r in records
        cid = Int(r["contiguous_index"]); zt = String(r["zone_type"])
        get!(id_to_types, cid, Set{String}()); push!(id_to_types[cid], zt)
    end
    zone_types_by_id = Dict(k => sort!(collect(v)) for (k,v) in id_to_types)
    ids_with_multiple = sort!([k for (k,v) in zone_types_by_id if length(v) > 1])
    return Dict(
        "num_unique_contiguous_indices"=>length(zone_types_by_id),
        "zone_types_by_id"=>zone_types_by_id,
        "ids_with_multiple_zone_types"=>ids_with_multiple,
    )
end

function export_to_glennht_conn(matches::Vector{Dict{String,Any}},
                                outer_faces::Vector{Dict{String,Int}},
                                filename::AbstractString,
                                gif_pairs::Vector{Dict{String,Int}},
                                gif_faces::Vector{Dict{String,Int}},
                                volume_zones::Vector{Dict{String,Any}})
    lines = IOBuffer()
    # matches
    nMatches = length(matches)
    @printf(lines, "%d\n", nMatches)
    for match in matches
        for side in ("block1","block2")
            blk = match[side]::Dict{String,Any}
            block_indx = Int(blk["block_index"])+1
            IMIN = Int(blk["IMIN"])+1; JMIN = Int(blk["JMIN"])+1; KMIN = Int(blk["KMIN"])+1
            IMAX = Int(blk["IMAX"])+1; JMAX = Int(blk["JMAX"])+1; KMAX = Int(blk["KMAX"])+1
            @printf(lines, "%3d\t%5d %5d %5d\t%5d %5d %5d\n", block_indx, IMIN,JMIN,KMIN, IMAX,JMAX,KMAX)
        end
    end

    # surfaces (outer + gif_faces), sorted by id
    outer = copy(outer_faces)
    append!(outer, gif_faces)
    sort!(outer; by=x->Int(get(x,"id",0)))
    @printf(lines, "%d\n", length(outer))
    for s in outer
        block_indx = Int(s["block_index"])+1
        IMIN = Int(s["IMIN"])+1; JMIN = Int(s["JMIN"])+1; KMIN = Int(s["KMIN"])+1
        IMAX = Int(s["IMAX"])+1; JMAX = Int(s["JMAX"])+1; KMAX = Int(s["KMAX"])+1
        sid  = Int(get(s,"id",0))
        @printf(lines, "%3d\t%5d %5d %5d\t%5d %5d %5d\t%4d\n", block_indx, IMIN,JMIN,KMIN, IMAX,JMAX,KMAX, sid)
    end

    # GIF pairs
    @printf(lines, "%d\n", length(gif_pairs))
    for p in gif_pairs
        a = Int(get(p,"a",0)); b = Int(get(p,"b",0))
        @printf(lines, "%d %d -2 1\n", a, b)
    end

    # volume zones
    summary = summarize_contiguous(volume_zones)
    @printf(lines, "%d\n", Int(summary["num_unique_contiguous_indices"]))
    # list zone ids on one line
    for (cid, _) in sort!(collect(summary["zone_types_by_id"]); by=first)
        @printf(lines, "%d ", cid)
    end
    write(lines, "\n")
    # print contiguous_id per block (10 per line, matching Python intent)
    cols = 10; col = 0
    for v in volume_zones
        cid = Int(v["contiguous_index"])
        col += 1
        if col == cols
            @printf(lines, "%d\n", cid); col = 0
        else
            @printf(lines, "%d ", cid)
        end
    end
    col != 0 && write(lines, "\n")

    out = ensure_extension(filename, ".ght_conn")
    open(out, "w") do io
        write(io, String(take!(lines)))
    end
    return nothing
end

