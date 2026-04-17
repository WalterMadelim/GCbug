"""
General master-subproblem that can be used to generate PI/SB/B-cut
"""
module General
import ..Case2383, ..Uvx, ..Settings, JuMP, Gurobi

_1(sub,s,m,T,t,a...) = setindex!(sub, sub!(s, T, t, m, a...), s)
_spawn_1(a...) = Threads.@spawn(_1(a...))
function model!(sub, t, T, tks, a...)
    inn = similar(tks, JuMP.Model)
    Settings.Model!(inn, tks)
    mst = mst!(t, Settings.Model(), a...)
    for (s,m)=enumerate(inn)
        tks[s] = _spawn_1(sub,s,m,T,t,a...)
    end
    foreach(wait, tks)
    mst
end

const NtType=@NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, qy::JuMP.AffExpr, Ci::Vector{Int32}, Cd::Vector{Float64}, Xl::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}, xllen::Int64, blen::Int64}
function sub!(s, T, t, m, LineˈP, F, WD, LD, GD, uH, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
    ReserveK01v, LoadKref, times,Times, lines = Case2383.Reserve_Curve[reserve_type], Case2383.Load_Curve, t+1:t+T,t:t+T, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    
    Common, Qˈy = _0(), _0()
    egp         =  Dict(i    =>_0() for i=Times)
    rgp         = [Dict(i    =>_0() for i=Times) for z=1:4] # 1≤z≤4
    pfe         =  Dict((l,i)=>_0() for i=Times  for l=lines)
    
    JuMP.@variables(m, begin
        # A minimal set of Linking variables from 1st-stage
        0 ≤ aL[z=eachindex(GDˈn), g=eachindex(GDˈn[z]); GDˈMmRatio[z][g] > 1.1] # (pA1_Link) to be fixed
        0 ≤ u1[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] # Link, to be fixed
        0 ≤ x1[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] ≤ 1 # (Int) Link, to be fixed
        # The rest of (non-Linking) Binary variables
        0 ≤ x2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=times] ≤ 1 # At first, relax Int-constr
        # The rest (less vital) variables, intra-1st-stage or intra-2nd-stage
        0 ≤ a1[z=eachindex(GDˈn), g=eachindex(GDˈn[z]); GDˈMmRatio[z][g] ≤ 1.1] # (pA1_nonLink)
        0 ≤ a2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=times] # ✅ p := (Pmin)x + pA
        0 ≤ u2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=times]
        0           ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), Times]
        0           ≤  ϖ[z=eachindex(WDˈN),                       Times] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn),                       Times] ≤ 1
        -LineˈP[l]  ≤ pf[l=lines,                                 Times] ≤ LineˈP[l]
    end)

    #= Wind/Load =# for i=Times
        egpˈi, ReserveK01, Obj = egp[i], ReserveK01v[i], ifelse(i===t, Common, Qˈy)
        for (z,Ø)=enumerate(WDˈN)
            Kˈz, ϖˈz, ϖCostˈz = windMatˈs[i,z], ϖ[z,i], length(Ø)ϖCost
            JuMP.add_to_expression!(Obj, ϖCostˈz, ϖˈz)
            for (g,node)=enumerate(Ø)
                pWind = Kˈz * @inbounds(WDˈPMax[z][g])
                JuMP.add_to_expression!(egpˈi, pWind)
                JuMP.add_to_expression!(egpˈi, -pWind, ϖˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l,i], @inbounds F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * pWind)
                    JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
                end
            end
        end
        for (z,Ø)=enumerate(LDˈn)
            ζCostˈz, ζˈz = length(Ø)ζCost, ζ[z,i] # load shed
            JuMP.add_to_expression!(Obj, ζCostˈz, ζˈz)
            1≤z≤4 && (rgpˈz = rgp[z][i])
            for (g,node)=enumerate(Ø)
                @inbounds begin
                    Ty, ConstLoad = LDˈi[z][g], LDˈg[z][g]
                    load = LoadKref[Ty][i] * ConstLoad
                end
                if 1≤z≤4
                    upResDemand = *(ReserveK01, ReserveK2D, ConstLoad)
                    JuMP.add_to_expression!(rgpˈz, -upResDemand)
                    JuMP.add_to_expression!(rgpˈz, upResDemand, ζˈz)
                end
                JuMP.add_to_expression!(egpˈi, -load)
                JuMP.add_to_expression!(egpˈi, load, ζˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l,i], @inbounds F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * -load)
                    JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
                end
            end
        end
    end

    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g]
            #= special references =# pAHˈzg = pAH[z][g] # ::Float64
            uHˈzg, xHˈzg = uH[z][g], xH[z][g] # ::AbstractVector
            u1ˈzg, x1ˈzg =  u1[z,g], x1[z,g] # ::JuMP.VariableRef
            pA1ˈzg = MmRatio > 1.1 ? aL[z,g] : a1[z,g] # ::JuMP.VariableRef
            rui, rdi, r2s = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            suC = GDˈsci[z][g] * linC
        end
        Gwid = PMax-Pmin
        Ru, Rd = Gwid/rui, Gwid/rdi
        Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
        # if UT == 1, use these coefficients
        T12 = Sd-PMax
        T13 = -max(Sd-Su,0.)
        T15 = -max(Su-Sd,0.)
        T16 = Su-PMax
        # Ramp constr, use these
        R1 = Su-Pmin
        R2 = Sd-Pmin-Rd
        let i=t, egpˈi=egp[i], urˈzg=ur[z,g,i], xHˈzgˈ1=xHˈzg[1]
            1≤z≤4 && JuMP.add_to_expression!(rgp[z][i], urˈzg)
            JuMP.add_to_expression!(egpˈi, Pmin, x1ˈzg)
            JuMP.add_to_expression!(egpˈi, pA1ˈzg)
            for l=lines
                (pfl, Cnl) = (pfe[l,i], @inbounds F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, x1ˈzg)
                JuMP.add_to_expression!(pfl, Cnl,       pA1ˈzg)
            end
            JuMP.add_to_expression!(Common, suC, u1ˈzg)
            JuMP.add_to_expression!(Common, linC * Pmin, x1ˈzg)
            JuMP.add_to_expression!(Common, linC, pA1ˈzg)
            #= U-V-X Polytope Constr =# JuMP.@constraints(m, begin
                u1ˈzg - x1ˈzg + xHˈzgˈ1 ≥ 0 # write this directly
                Uvx.u!(-x1ˈzg,UT,i,t,z,g,uHˈzg,u1ˈzg,u2) ≤ 0
                Uvx.v!(_0(),DT,i,t,z,g,xHˈzg,x1ˈzg,x2,uHˈzg,u1ˈzg,u2) ≤ 1
            end)
            if MmRatio > 1.1
                JuMP.@constraints(m, begin
                    urˈzg+pA1ˈzg          ≤ Gwid*x1ˈzg + T16*u1ˈzg
                    urˈzg+pA1ˈzg - pAHˈzg ≤ Ru * xHˈzgˈ1 + R1 *  u1ˈzg
                    pAHˈzg      -  pA1ˈzg ≤ Rd * xHˈzgˈ1 + R2 * (u1ˈzg - x1ˈzg + xHˈzgˈ1)
                end)
            else # add a naive ub-output limit, but no ramp limits
                JuMP.@constraint(m, urˈzg+pA1ˈzg ≤ Gwid*x1ˈzg) # This is THE naive
            end
        end
        for i=times # 2nd-Linking-to-1st _or_ intra-2nd-stage
            urˈzg, egpˈi, xˈzg, uˈzg, pAˈzg = ur[z,g,i], egp[i], x2[z,g,i], u2[z,g,i], a2[z,g,i]
            1≤z≤4 && JuMP.add_to_expression!(rgp[z][i], urˈzg) 
            JuMP.add_to_expression!(egpˈi, Pmin, xˈzg)
            JuMP.add_to_expression!(egpˈi, pAˈzg)
            for l=lines
                (pfl, Cnl) = (pfe[l,i], @inbounds F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, xˈzg)
                JuMP.add_to_expression!(pfl, Cnl,       pAˈzg)
            end
            JuMP.add_to_expression!(Qˈy, suC, uˈzg)
            JuMP.add_to_expression!(Qˈy, linC * Pmin, xˈzg)
            JuMP.add_to_expression!(Qˈy, linC, pAˈzg)
            #= U-V-X Polytope Constr =# JuMP.@constraints(m, begin
                Uvx.vge0!(1.0*uˈzg,i,t,z,g,x1ˈzg,x2) ≥ 0 # ✅ 0 ≤ v[1] := u[1]-x[1]+x[0]
                Uvx.u!(-xˈzg,UT,i,t,z,g,uHˈzg,u1ˈzg,u2) ≤ 0
                Uvx.v!(_0(),DT,i,t,z,g,xHˈzg,x1ˈzg,x2,uHˈzg,u1ˈzg,u2) ≤ 1
            end)
            if MmRatio > 1.1
                #= UpperLimit Constr =# if i == t+T
                    JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg + T16 * uˈzg) # This is THE classic
                elseif UT > 1
                    JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg + T12 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T16 * uˈzg)
                else
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg ≤ Gwid*xˈzg + T12 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T13 * uˈzg
                        urˈzg+pAˈzg ≤ Gwid*xˈzg + T15 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T16 * uˈzg
                    end)
                end
                #= Normal Ramp Constr =# if i == t+1
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg - pA1ˈzg ≤ Ru * x1ˈzg + R1 *  uˈzg
                        pA1ˈzg      - pAˈzg  ≤ Rd * x1ˈzg + R2 * (uˈzg - xˈzg + x1ˈzg)
                    end)
                else
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg - a2[z,g,i-1] ≤ Ru * x2[z,g,i-1] + R1 *  uˈzg
                        a2[z,g,i-1] - pAˈzg       ≤ Rd * x2[z,g,i-1] + R2 * (uˈzg - xˈzg + x2[z,g,i-1])
                    end)
                end
            else # add a naive ub-output limit, but no ramp limits
                JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg) # This is THE naive
            end
        end
    end

    JuMP.@objective(m, Min, Common + Qˈy) # For a certain scene, assign an initial objfun
    JuMP.@constraints(m, begin
        [i=Times],          egp[i] == 0
        [z=1:4, i=Times],   rgp[z][i] ≥ 0
        [l=lines,i=Times],  pfe[l,i] == pf[l,i]
    end)

    #=Also bst=# xllen = length(aL)+length(u1)+length(x1)
    blen = length(x2)
    Ci = Cint[range(length(WDˈS)+1; length=xllen); s]
    Cd = similar(Ci, Cdouble); Cd[end] = -1. # Π'x - θs ≤ -Cn || -Π'x + θs ≥ Cn
    Xl = view(Cd, 1:xllen)
    o = m.moi_backend
    refd = Ref{Cdouble}()
    refi = Ref{Cint}()
    qy = Qˈy
    (; m, o, refd, refi, qy, Ci, Cd, Xl, xllen, blen)
end

function mst!(t, m, LineˈP, F, WD, LD, GD, uH, xH, pAH; s=1, ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
    ReserveK01, LoadKref, lines = Case2383.Reserve_Curve[reserve_type][t], Case2383.Load_Curve, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    
    Common      = _0()
    egp         = _0()
    rgp         = [_0() for z=1:4] # 1≤z≤4
    pfe         = Dict(l=>_0() for l=lines)

    JuMP.@variables(m, begin
        common # This _variable_ can be used to construct PI-cuts
        θ[eachindex(WDˈS)]
        # A minimal set of Linking variables from 1st-stage
        0 ≤ aL[z=eachindex(GDˈn), g=eachindex(GDˈn[z]); GDˈMmRatio[z][g] > 1.1] # (pA1_Link) to be fixed
        0 ≤ u1[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] # Link, to be fixed
        0 ≤ x1[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] ≤ 1 # (Int) Link, to be fixed
        # The rest (less vital) variables, intra-1st-stage or intra-2nd-stage
        0 ≤ a1[z=eachindex(GDˈn), g=eachindex(GDˈn[z]); GDˈMmRatio[z][g] ≤ 1.1] # (pA1_nonLink)
        0           ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z])]
        0           ≤  ϖ[z=eachindex(WDˈN)                      ] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn)                      ] ≤ 1
        -LineˈP[l]  ≤ pf[l=lines                                ] ≤ LineˈP[l]
    end);

    #= Wind/Load =# let i=t
        for (z,Ø)=enumerate(WDˈN)
            Kˈz, ϖˈz, ϖCostˈz = windMatˈs[i,z], ϖ[z], length(Ø)ϖCost
            JuMP.add_to_expression!(Common, ϖCostˈz, ϖˈz)
            for (g,node)=enumerate(Ø)
                pWind = Kˈz * @inbounds(WDˈPMax[z][g])
                JuMP.add_to_expression!(egp, pWind)
                JuMP.add_to_expression!(egp, -pWind, ϖˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l], @inbounds F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * pWind)
                    JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
                end
            end
        end
        for (z,Ø)=enumerate(LDˈn)
            ζCostˈz, ζˈz = length(Ø)ζCost, ζ[z] # load shed
            JuMP.add_to_expression!(Common, ζCostˈz, ζˈz)
            1≤z≤4 && (rgpˈz = rgp[z])
            for (g,node)=enumerate(Ø)
                @inbounds begin
                    Ty, ConstLoad = LDˈi[z][g], LDˈg[z][g]
                    load = LoadKref[Ty][i] * ConstLoad
                end
                if 1≤z≤4
                    upResDemand = *(ReserveK01, ReserveK2D, ConstLoad)
                    JuMP.add_to_expression!(rgpˈz, -upResDemand)
                    JuMP.add_to_expression!(rgpˈz, upResDemand, ζˈz)
                end
                JuMP.add_to_expression!(egp, -load)
                JuMP.add_to_expression!(egp, load, ζˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l], @inbounds F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * -load)
                    JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
                end
            end
        end
    end

    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g]
            #= special references =# pAHˈzg = pAH[z][g] # ::Float64
            uHˈzg, xHˈzg = uH[z][g], xH[z][g] # ::AbstractVector
            u1ˈzg, x1ˈzg =  u1[z,g], x1[z,g] # ::JuMP.VariableRef
            pA1ˈzg = MmRatio > 1.1 ? aL[z,g] : a1[z,g] # ::JuMP.VariableRef
            rui, rdi, r2s = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            suC = GDˈsci[z][g] * linC
        end
        Gwid = PMax-Pmin
        Ru, Rd = Gwid/rui, Gwid/rdi
        Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
        # if UT == 1, use these coefficients
        T12 = Sd-PMax
        T13 = -max(Sd-Su,0.)
        T15 = -max(Su-Sd,0.)
        T16 = Su-PMax
        # Ramp constr, use these
        R1 = Su-Pmin
        R2 = Sd-Pmin-Rd
        let i=t, urˈzg=ur[z,g], xHˈzgˈ1=xHˈzg[1], u2=nothing, x2=nothing
            1≤z≤4 && JuMP.add_to_expression!(rgp[z], urˈzg)
            JuMP.add_to_expression!(egp, Pmin, x1ˈzg)
            JuMP.add_to_expression!(egp, pA1ˈzg)
            for l=lines
                (pfl, Cnl) = (pfe[l], @inbounds F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, x1ˈzg)
                JuMP.add_to_expression!(pfl, Cnl,       pA1ˈzg)
            end
            JuMP.add_to_expression!(Common, suC, u1ˈzg)
            JuMP.add_to_expression!(Common, linC * Pmin, x1ˈzg)
            JuMP.add_to_expression!(Common, linC, pA1ˈzg)
            #= U-V-X Polytope Constr =# JuMP.@constraints(m, begin
                u1ˈzg - x1ˈzg + xHˈzgˈ1 ≥ 0 # write this directly
                Uvx.u!(-x1ˈzg,UT,i,t,z,g,uHˈzg,u1ˈzg,u2) ≤ 0
                Uvx.v!(_0(),DT,i,t,z,g,xHˈzg,x1ˈzg,x2,uHˈzg,u1ˈzg,u2) ≤ 1
            end)
            if MmRatio > 1.1
                JuMP.@constraints(m, begin
                    urˈzg+pA1ˈzg          ≤ Gwid*x1ˈzg + T16*u1ˈzg
                    urˈzg+pA1ˈzg - pAHˈzg ≤ Ru * xHˈzgˈ1 + R1 *  u1ˈzg
                    pAHˈzg      -  pA1ˈzg ≤ Rd * xHˈzgˈ1 + R2 * (u1ˈzg - x1ˈzg + xHˈzgˈ1)
                end)
            else # add a naive ub-output limit, but no ramp limits
                JuMP.@constraint(m, urˈzg+pA1ˈzg ≤ Gwid*x1ˈzg) # This is THE naive
            end
        end
    end

    JuMP.@objective(m, Min, Common + sum(θ)/length(θ)) # It is set only once
    JuMP.@constraints(m, begin
        Common == common
        egp == 0
        rgp .≥ 0
        [l=lines],  pfe[l] == pf[l]
    end)
    xlst = 1+length(θ)
    xllen = length(aL)+length(u1)+length(x1)
    Xl = Vector{Cdouble}(undef, xllen)
    o = m.moi_backend
    refd = Ref{Cdouble}()
    refi = Ref{Cint}()
    bst, blen = sum(length, (common,θ,aL,u1)), length(x1)
    Θ = similar(θ, Cdouble)
    (; m, o, refd, refi, Θ, xlst, xllen, Xl, bst, blen)
end

_0() = JuMP.AffExpr(0.)
end

