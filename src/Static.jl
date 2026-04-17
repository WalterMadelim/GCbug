module Static
import JuMP, ..Case2383, ..Settings

_0() = JuMP.AffExpr(0.)
function Model2383!(m, t, F, liD, WD, LD, GN, GMax, Gmin;ReserveK2D=0.04, reserve_type=1)
    ReserveK01, LoadKref = Case2383.Reserve_Curve[reserve_type][t], Case2383.Load_Curve
    LDˈi, LDˈg, LDˈn = LD.i, LD.g, LD.n
    lines = keys(liD)
    pfe, egp, rgp = Dict(k=>_0() for k=lines), _0(), [_0() for z=1:4]
    JuMP.@variables(m, begin
        -liD[k] ≤ pf[k=lines] ≤ liD[k]
        x[z=eachindex(GN), g=eachindex(GN[z])], Bin # static problem, no turn-on/off actions
        0 ≤ pA[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ur[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ϖ[z=eachindex(WD.N)] ≤ 0
        0 ≤ ζ[z=eachindex(LDˈn)] ≤ 0
    end); JuMP.@constraint(m, [z=eachindex(GN), g=eachindex(GN[z])], 
        pA[z,g]+ur[z,g] ≤ (GMax[z][g]-Gmin[z][g])x[z,g])
    for (z,v)=enumerate(LDˈg)
        ζˈz = ζ[z] # load shed
        1≤z≤4 && (rˈz = rgp[z]) 
        for (g,p)=enumerate(v)
            Ty, n = LDˈi[z][g], LDˈn[z][g]
            load = LoadKref[Ty][t] * p
            JuMP.add_to_expression!(egp, -load)
            JuMP.add_to_expression!(egp, load, ζˈz)
            if 1≤z≤4
                upResDemand = ReserveK01 * ReserveK2D * p
                JuMP.add_to_expression!(rˈz, -upResDemand)
                JuMP.add_to_expression!(rˈz, upResDemand, ζˈz)
            end
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * -load)
                JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
            end
        end
    end
    for (z,v)=enumerate(WD.PMax)
        Kˈz, ϖˈz = WD.S[1][1, z], ϖ[z]
        for (g,pMa)=enumerate(v)
            pWind = Kˈz * pMa
            JuMP.add_to_expression!(egp, pWind)
            JuMP.add_to_expression!(egp, -pWind, ϖˈz)
            n = WD.N[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pWind)
                JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
            end
        end
    end
    for (z,v)=enumerate(Gmin)
        1≤z≤4 && (rˈz = rgp[z])
        for (g,pmin)=enumerate(v)
            1≤z≤4 && JuMP.add_to_expression!(rˈz, ur[z,g])
            JuMP.add_to_expression!(egp, pmin, x[z,g])
            JuMP.add_to_expression!(egp, pA[z,g])
            n = GN[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pmin, x[z,g])
                JuMP.add_to_expression!(pfl, Cnl, pA[z,g])
            end
        end
    end
    JuMP.@constraints(m, begin
        egp == 0
        rgp .>= 0
        [k=lines], pfe[k] == pf[k]
    end)
    JuMP.@objective(m, Min, sum(rand(5.0:1e-4:8.0)i for i=x) + sum(rand(0.5:1e-4:1.5)i for i=pA))
end
function get_History_2383(t, F, Line, WD, LD, GD)
    GDˈUT = GD.UT
    uL = max(maximum(maximum, GDˈUT),maximum(maximum, GD.DT))-1 # uH is not only for `u` but also for `v`
    uH =  [[falses(uL  ) for _=v] for v=GDˈUT]
    xH =  [[falses(uL+1) for _=v] for v=GDˈUT]
    pAH = [[NaN          for _=v] for v=GDˈUT]
    model = Settings.Model()
    Model2383!(model, t, F, Line.P, WD, LD, GD.n, GD.pmax, GD.pmin)
    m = (m = model, o = model.moi_backend, refd = Ref{Cdouble}(), refi = Ref{Cint}())
    ter = Settings.opt_and_ter(m)
    ter == 2 || error("terminate code $ter")
    pA, x = m.m[:pA], m.m[:x]
    for (z,v)=enumerate(uH), (g,bv)=enumerate(v)
        pAH[z][g] = Settings.getxdblattrelement(m, pA[z,g], "X")
        if round(Bool, Settings.getxdblattrelement(m, x[z,g], "X"))
            if rand() ≤ .9
                j = rand(eachindex(bv))
                xH[z][g][1:j] .= bv[j] = true
            else
                xH[z][g] .= true
            end
        else
            if rand() ≤ .9
                j = rand(eachindex(bv))
                xH[z][g][j+1:end] .= true
            end
        end
    end
    xH,uH,pAH
end

# Below are not case2383
function Model!(m, s, i, F, GD, WD, DD; PfMax, KupRe) # This is a simplified single-period model
    GDˈi, GDˈD, GDˈn = GD.i, GD.D, GD.n
    wMat, WDˈN, WDˈPMax, WDˈZone = WD.D[s], WD.N, WD.PMax, WD.Zone
    DDˈi, DDˈD, DDˈn, DDˈC = DD.i, DD.D, DD.n, DD.C
    lines = 1:size(F,1)
    rgp, egp, fgp = _0(), _0(), [_0() for l=lines] # sys-level, per time slot
    JuMP.@variables(m, begin
        -PfMax ≤ pf[l=lines] ≤ PfMax
        0. ≤ ϖ[z=WDˈZone, g=eachindex(WDˈN[z])] ≤ WDˈPMax[z][g]wMat[i,z]
        x[GDˈi], Bin
        0. ≤ pA[GDˈi]
        0. ≤ ur[GDˈi]
        p[GDˈi]
    end)
    for (type, D0, node) = zip(DDˈi, DDˈD, DDˈn)
        Dmd = DDˈC[type][i]D0
        JuMP.add_to_expression!(rgp,  KupRe * -Dmd)
        JuMP.add_to_expression!(egp, -Dmd)
        for l=lines JuMP.add_to_expression!(fgp[l], F[l,node] * -Dmd) end
    end
    for z=WDˈZone, (g, node)=enumerate(WDˈN[z])
        ϖUB, ϖzg = WDˈPMax[z][g]wMat[i,z], ϖ[z,g]
        JuMP.add_to_expression!(egp, ϖUB)
        JuMP.add_to_expression!(egp, -1., ϖzg)
        for l=lines
            el, Cnl = fgp[l], F[l,node]
            JuMP.add_to_expression!(el, Cnl * ϖUB)
            JuMP.add_to_expression!(el, -Cnl, ϖzg)
        end
    end
    for (g, GD)=GDˈD
        node = GDˈn[g]
        JuMP.add_to_expression!(rgp, ur[g])
        JuMP.add_to_expression!(egp, p[g])
        for l=lines
            el, Cnl = fgp[l], F[l,node]
            JuMP.add_to_expression!(el, Cnl, p[g])
        end
        # local constrs
        (GMin, GMax) = GD.o; GWid = GMax-GMin
        JuMP.@constraint(m, ur[g]+pA[g] ≤ (GWid)x[g])
        JuMP.@constraint(m, p[g] == (GMin)x[g] + pA[g])
    end
    JuMP.@constraints(m, begin # sys-level, must goes after `add_to_expression!`
        rgp >= 0.
        egp == 0.
        [l=lines], fgp[l] == pf[l]
    end)
end

function getHistory(s,t,F,GD,WD,DD,uvHistory;PfMax=6.47,KupRe=0.1)
    m = Settings.Model()
    Model!(m, s, t, F, GD, WD, DD; PfMax=PfMax, KupRe=KupRe)
    JuMP.optimize!(m)
    JuMP.termination_status(m) === JuMP.OPTIMAL || error()
    x,pA = m[:x],m[:pA]
    for g=GD.i
        if JuMP.value(x[g]) < 0.5
            rand()>0.1 && (uvHistory[g, 1, rand(Gen.HisRg(t))]=true)
        else
            rand()>0.1 && (uvHistory[g, 2, rand(Gen.HisRg(t))]=true)
        end
    end
    xHistory, pAHistory = JuMP.value.(x), JuMP.value.(pA)
end

end


