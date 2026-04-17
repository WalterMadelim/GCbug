"""
helper functions for adding the UVX-Polytope
    0 ≤ v[1] := u[1]-x[1]+x[0]
    &&
    sum(u) ≤ x
    &&
    sum(v) <= 1     - x
    v[i]   = u[i]   - x[i]   + x[i-1]
    v[i-1] = u[i-1] - x[i-1] + x[i-2]
Conclusion: you need to add a series of `u`-counterpart
Finally there's a `x` with lag one
"""
module Uvx
import JuMP

pickAdd!(e,i,t,z,g,uH::AbstractVector,u1::JuMP.VariableRef,u) = if i > t
    JuMP.add_to_expression!(e, u[z,g,i])
elseif t > i
    JuMP.add_to_expression!(e, uH[t-i])
else
    JuMP.add_to_expression!(e, u1)
end
function u!(e, Len, i, a...) # sum(u) ≤ x
    cnt = 0
    while true # the initial `i` is at an outer layer
        cnt == Len && return e
        pickAdd!(e, i, a...)
        i -= 1 # here we reuse `i`, which is a different concept
        cnt += 1
    end
end
function v!(e, Len, i, t, z,g, xH,x1,x, uH,u1,u)
    cnt = 0
    while true # the initial `i` is at an outer layer
        cnt == Len && break
        pickAdd!(e,i,t,z,g,uH,u1,u)
        i -= 1 # here we reuse `i`, which is a different concept
        cnt += 1
    end
    pickAdd!(e,i,t,z,g,xH,x1,x)
    e
end
function vge0!(e,i,t,z,g,x1,x) # 0 ≤ v[1] := u[1]-x[1]+x[0]
    JuMP.add_to_expression!(e, -1., x[z,g,i])
    i -= 1
    if i == t
        JuMP.add_to_expression!(e, x1)
    elseif i > t
        JuMP.add_to_expression!(e, x[z,g,i])
    else
        error()
    end
    e
end

end
