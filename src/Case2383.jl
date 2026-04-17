module Case2383
import PowerModels, PGLib
const CaseName = "pglib_opf_case2383wp_k.m"
const Reserve_Curve=(
    [0.72, 0.71, 0.72, 0.73, 0.76, 0.84, 0.91, 0.91, 0.91, 0.9, 0.9, 0.89, 0.89, 0.88, 0.87, 0.86, 0.88, 0.97, 1.0, 0.98, 0.94, 0.88, 0.81, 0.75],
    [0.74, 0.73, 0.73, 0.73, 0.77, 0.84, 0.9, 0.9, 0.89, 0.9, 0.9, 0.9, 0.9, 0.89, 0.87, 0.87, 0.88, 0.97, 1.0, 0.98, 0.94, 0.88, 0.81, 0.77]
)
const Load_Curve = (
    [0.54, 0.53, 0.52, 0.52, 0.53, 0.57, 0.64, 0.68, 0.7, 0.71, 0.69, 0.67, 0.68, 0.71, 0.72, 0.73, 0.74, 0.77, 0.83, 0.81, 0.74, 0.68, 0.62, 0.57],
    [0.45, 0.44, 0.44, 0.44, 0.44, 0.49, 0.6, 0.69, 0.76, 0.8, 0.78, 0.78, 0.78, 0.79, 0.79, 0.75, 0.71, 0.65, 0.58, 0.53, 0.49, 0.47, 0.45, 0.44],
    [0.83, 0.83, 0.83, 0.83, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.9, 0.9, 0.9, 0.89, 0.88, 0.87, 0.86, 0.86, 0.85, 0.85, 0.84, 0.84],
    [0.5, 0.5, 0.5, 0.5, 0.51, 0.58, 0.7, 0.78, 0.8, 0.83, 0.79, 0.73, 0.72, 0.7, 0.69, 0.73, 0.78, 0.83, 0.75, 0.63, 0.58, 0.55, 0.53, 0.52],
    [0.59, 0.56, 0.54, 0.53, 0.54, 0.59, 0.64, 0.66, 0.69, 0.71, 0.74, 0.76, 0.79, 0.81, 0.84, 0.81, 0.76, 0.69, 0.64, 0.61, 0.59, 0.64, 0.71, 0.79],
    [0.53, 0.52, 0.52, 0.52, 0.52, 0.58, 0.68, 0.76, 0.81, 0.83, 0.82, 0.81, 0.81, 0.81, 0.8, 0.77, 0.72, 0.68, 0.63, 0.61, 0.58, 0.57, 0.56, 0.54]
)

get_Case_Dict()=(PowerModels.silence();PowerModels.make_basic_network(PGLib.pglib(CaseName)))
_0!(F; a=9e-7)=for (i,e)=enumerate(F) if -a<e<a setindex!(F, 0., i) end end
get_PTDF(CaD)=(F=PowerModels.calc_basic_ptdf_matrix(CaD);_0!(F);F)

struct Bbyz
    v::Vector{Vector{Int}} # zone relations
    m::Matrix{Vector{Int}} # branch indicies, organized into zone pair
    P::Dict{Int, Float64} # explicitly extract PfMax for type stability
    Bbyz(N)=new([Int[] for _=1:N-1], Matrix{Vector{Int}}(undef,N,N), Dict{Int, Float64}())
end
function branch_by_zone(CaD; PfK=1.0)
    o = Bbyz(6)
    CaDˈbus, CaDˈbr = CaD["bus"], CaD["branch"]
    for i=1:length(CaDˈbr)
        CaDˈbrˈi = CaDˈbr[string(i)]
        f,t = CaDˈbrˈi["f_bus"],CaDˈbrˈi["t_bus"]
        f,t = CaDˈbus[string(f)]["zone"],CaDˈbus[string(t)]["zone"]
        f == t && continue
        o.P[i] = PfK * CaDˈbrˈi["rate_a"] # This original Data is usable. To tighten it, setting the multiplier to 0.6 would be too tight to be feasible.
        f,t = min(f,t), max(f,t)
        if t in o.v[f]
            push!(o.m[f,t], i)
        else
            o.m[f,t], _ = [i], push!(o.v[f], t)
        end
    end
    foreach(sort!, o.v)
    o
end
function get_load(CaD)
    N, G, I = [Int[] for _=1:6], [Float64[] for _=1:6], [Int[] for _=1:6] # `I` is index of Curve
    CaDˈbus,CaDˈload = CaD["bus"],CaD["load"]
    Threshold = 0.001
    for g=1:length(CaDˈload)
        CaDˈloadˈg = CaDˈload[string(g)]
        p = CaDˈloadˈg["pd"]
        p < Threshold && continue
        n = CaDˈloadˈg["load_bus"]
        z = CaDˈbus[string(n)]["zone"]
        curve_index = rand(1:6) # It happens that #Zone == #LoadCurve == 6
        _,_,_ = push!(N[z], n),push!(G[z], p),push!(I[z], curve_index)
    end
    (n=N, g=G, i=I)
end

get_type(PMa) = if PMa < 0.5 # used by `case2383`
    i = 6 # gasCT (Smallest)
elseif PMa < 0.505
    i = 5 # oil
elseif PMa < 0.9
    i = 4 # hydro
elseif PMa < 1.4
    i = 3 # gasCC
elseif PMa < 10.
    i = 2 # coal
else
    i = 1 # nuclear (Biggest)
end

"""
Hold-time generator, same for Dn/Up Times
"""
_HT(Ty) = if Ty == 1 # nuclear
    23
elseif Ty == 2 # coal
    rand(12:23)
elseif Ty == 3 # gasCC
    rand(6:23)
elseif Ty == 4 # hydro
    rand(1:6)
elseif Ty == 5 # oil
    rand(3:12)
elseif Ty == 6 # gasCT
    rand(1:4)
end

"""
(Dn, UP)
This is Int, not readily (RD, RU). They depends on GMax, GMin
"""
_RampPeriods(Ty) = if Ty == 1 # nuclear
    23,23
elseif Ty == 2 # coal
    rand(6:17),rand(8:23)
elseif Ty == 3 # gasCC
    rand(2:8),rand(3:10)
elseif Ty == 4 # hydro
    2,2
elseif Ty == 5 # oil
    rand(2:6), rand(2:6)
elseif Ty == 6 # gasCT
    2,2
end

"""
relative_cost::Int
"""
_SUcostOprcostRatio(Ty) = if Ty == 1 # nuclear
    rand(50:150)
elseif Ty == 2 # coal
    rand(20:80)
elseif Ty == 3 # gasCC
    rand(10:40)
elseif Ty == 4 # hydro
    rand(1:5)
elseif Ty == 5 # oil
    rand(10:50)
elseif Ty == 6 # gasCT
    rand(2:15)
end

"""
SD/SU rate is calculated from RD/RU as clamp(R2S * Rrate, Gmin, GMax)
"""
_r2s() = rand(.8:.1:1.5)
_Ci() = [Int[] for _=1:6]
_Cf() = [Float64[] for _=1:6]
function get_gen(CaD)
    N, RUI, RDI, HUT, HDT, SCI = _Ci(), _Ci(), _Ci(), _Ci(), _Ci(), _Ci()
    G, MmRatio, Gmin, Clin, R2S = _Cf(), _Cf(), _Cf(), _Cf(), _Cf()
    CaDˈbus,CaDˈgen = CaD["bus"],CaD["gen"]
    Threshold = 0.001 * sum(v["pmax"] for (_,v)=CaDˈgen)
    for g=1:length(CaDˈgen)
        CaDˈgenˈg = CaDˈgen[string(g)]
        p = CaDˈgenˈg["pmax"]
        p < Threshold && continue
        n = CaDˈgenˈg["gen_bus"]
        z = CaDˈbus[string(n)]["zone"]
        push!(N[z], n) # node
        push!(G[z], p) # PMax
        Ty = get_type(p) # Type TODO need to export this info?
        pm = CaDˈgenˈg["pmin"]; push!(Gmin[z], pm) # Pmin
        ra = p / pm; push!(MmRatio[z], ra) # PMax/Pmin-Ratio 
        push!(HUT[z], _HT(Ty)); push!(HDT[z], _HT(Ty)) # minimun Dn/Up time
        C = CaDˈgenˈg["cost"][2]/5000; C < 0.2 && (C=rand(0.3:1e-4:0.5)); push!(Clin[z], C) # LinCost
        push!(SCI[z], _SUcostOprcostRatio(Ty)) # Start-Up Cost (indirect)
        ø1,ø2 = _RampPeriods(Ty); push!(RDI[z], ø1); push!(RUI[z], ø2) # normal Ramp (indirect)
        push!(R2S[z], _r2s()) # Start/Shut Ramp (indirect)
    end
    (
        n = N,
        pmin = Gmin,
        pmax = G,
        MmRatio = MmRatio, # if <1.1, don't add ramp constrs
        UT = HUT, # ::Int
        DT = HDT,
        C = Clin,
        sci = SCI, # start-up cost::Int
        rdi = RDI, # ramp-down::Int
        rui = RUI,
        r2s = R2S # a ratio::Float64
    )
end




end
