include(joinpath("src", "Settings.jl"));
include(joinpath("src", "Case2383.jl"));
include(joinpath("src", "WindGen.jl"));
include(joinpath("src", "Static.jl"));
include(joinpath("src", "Uvx.jl"));
include(joinpath("src", "General.jl"));

import JuMP, Gurobi, Random
Random.seed!(hash(1))
function build_2ssp!(sub, t, T, tks)
    WD = WindGen.get_case2383(tks, t)
    CaD = Case2383.get_Case_Dict()
    F = Case2383.get_PTDF(CaD)
    LD = Case2383.get_load(CaD)
    GD = Case2383.get_gen(CaD)
    Line = Case2383.branch_by_zone(CaD; PfK=0.5)
    xH,uH,pAH = Static.get_History_2383(t, F, Line, WD, LD, GD)
    println(" here, before >>")
    @time General.model!(sub, t, T, tks, Line.P, F, WD, LD, GD, uH, xH, pAH)
end;
function multi_th_opt(tks, sub)
    for (s, m)=enumerate(sub)
        tks[s] = Threads.@spawn(Settings.opt_ass_opt(m))
    end
    foreach(wait, tks)
end;
function pi_bound!(mst, sub, tks) # at the same time recover sub to BenCut state
    M = mst.m
    common, θ = M[:common], M[:θ] 
    for (s,m)=enumerate(sub)
        Cn = Settings.getmodeldblattr(m, "ObjBound")
        JuMP.@constraint(M, common + θ[s] ≥ Cn)
    end
    for (s,m)=enumerate(sub)
        tks[s] = Threads.@spawn(JuMP.set_objective_function(m.m, m.qy))
    end
    Settings.opt_ass_opt(mst)
    lb = Settings.getmodeldblattr(mst, "ObjBound")
    println("Perfect_Info_LB = $lb")
    foreach(wait, tks)
end;
function get_trial!(mst)
    ter = Settings.opt_and_ter(mst)
    if ter != 2
        JuMP.set_attribute(mst.m, "DualReductions", 0)
        JuMP.unset_silent(mst.m)
        ter = Settings.opt_and_ter(mst)
        error("mst Status code = $ter")
    end
    Gurobi.GRBgetdblattrarray(mst.o, "X", mst.xlst, mst.xllen, mst.Xl)
    Settings.getmodeldblattr(mst, "ObjBound"), Settings.getxdblattrelement(mst, 0, "X") # lb, Common
end;
function multi_solve(tks, sub, mst)
    Xl = mst.Xl
    for (s,n)=enumerate(sub)
        tks[s] = Threads.@spawn(scene!(n,Xl))
    end
    Gurobi.GRBgetdblattrarray(mst.o, "X", 1, mst.xllen, mst.Θ)
    foreach(wait, tks)
end;
function scene!(n,Xl)
    o, l = n.o, n.xllen
    Gurobi.GRBsetdblattrarray(o, "LB", 0, l, Xl)
    Gurobi.GRBsetdblattrarray(o, "UB", 0, l, Xl)
    Settings.opt_and_ter(n) == 2 || error()
end;
function bwd_in_sequential(sub, mst, ub)
    vioCnt, vioSum, S = 0, 0., length(sub)
    for (n, Θs)=zip(sub, mst.Θ)
        success, vio, obj = _scene_bwd(n, Θs, mst)
        ub += obj/S
        success || continue
        vioCnt += 1
        vioSum += vio
    end
    vioCnt === 0 && return (false, 0., ub) # no cut gened
    (true, vioSum/vioCnt, ub) # needs further train
end;
function _scene_bwd(n, Θs, mst)
    obj = Settings.getmodeldblattr(n, "ObjVal")
    vio = obj - Θs
    success = vio > 1e-4
    if success
        l, Cd = n.xllen, n.Cd
        Gurobi.GRBgetdblattrarray(n.o, "RC", 0, l, Cd)
        Gn = n.Xl'mst.Xl - obj # pure Benders' formula
        Gurobi.GRBaddconstr(mst.o, l+1, n.Ci, Cd, Cchar('<'), Gn, C_NULL)
    end
    success, vio, obj
end;

S = 3; # 256*20 already took >428GB
t,T = 1,23;
tks = Vector{Task}(undef, S);
sub = similar(tks, General.NtType);
mst = build_2ssp!(sub, t, T, tks);

multi_th_opt(tks, sub)
pi_bound!(mst, sub, tks)
lb, common = get_trial!(mst)

for k=1:20
    multi_solve(tks, sub, mst)
    proceed, vioMean, ub = bwd_in_sequential(sub, mst, common)
    lb, common = get_trial!(mst)
    agap = ub - lb
    rgap = agap/ub
    println("lb = $lb, agap = $agap, vioMean=$vioMean, rgap = $rgap")
end

tks = Vector{Task}(undef, S);
sub = similar(tks, General.NtType);
mst = build_2ssp!(sub, t, T, tks);
