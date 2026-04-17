module Settings
import JuMP, Gurobi
_gcc(o,x) = Gurobi.c_column(o,JuMP.index(x))
const C = Dict{String, Any}("OutputFlag"=>0,"Threads"=>1) # "SoftMemLimit" => 4e2

# m::MyNamedTuple
# |-m::Model
# |-o::Optimizer
# |-refd::ref_Cdouble
# |-refi::ref_Cint

# ✅ e.g. setting `x::JuMP.VariableRef`'s "Obj" Attribute
getxdblattrelement(m, i::Integer, str) = (r=m.refd;Gurobi.GRBgetdblattrelement(m.o, str,i,r);r.x)
getxdblattrelement(m, x::JuMP.VariableRef, str) = getxdblattrelement(m, _gcc(m.o, x), str)
setxdblattrelement(m, i::Integer, str, v) = Gurobi.GRBsetdblattrelement(m.o, str, i, v)
setxdblattrelement(m, x::JuMP.VariableRef, str, v) = setxdblattrelement(m, _gcc(m.o, x), str, v)

# ✅ e.g. add a constr
# Gurobi.GRBaddconstr(o, read_len, ci, cd, Cchar('=' or '>'), Gn, C_NULL)

# ✅ e.g. query "X" of the trial vec in 1st-stage, and then fix the copy vec in 2nd-stage
# Gurobi.GRBgetdblattrarray(o, "X", start, len, Ptr)
# Gurobi.GRBsetdblattrarray(o, "LB", start, len, Ptr)
# Gurobi.GRBsetcharattrarray(o, "VType", start, len, fill(Cchar('B'), len)) # ✅ e.g. set variables to binary

# ✅ e.g. query ObjVal
function getmodeldblattr(m, str)
    r = m.refd
    Gurobi.GRBgetdblattr(m.o, str, r)
    r.x
end
# ✅ e.g. query NumConstrs
function getmodelintattr(m, str)
    r = m.refi
    Gurobi.GRBgetintattr(m.o, str, r)
    r.x
end

opt_ass_opt(m) = (opt_and_ter(m)==2 || error())
function opt_and_ter(m)
    o, r = m.o, m.refi
    Gurobi.GRBoptimize(o)
    Gurobi.GRBgetintattr(o, "Status", r)
    r.x
end

Env() = Gurobi.Env(C) # generate a _new_ one as defined by `Gurobi.Env`
function Model()
    m = JuMP.direct_model(Gurobi.Optimizer(Env()))
    JuMP.set_string_names_on_creation(m, false)
    m
end
function Model!(v, t #=tks=#)
    for (j,k)=zip(eachindex(v),eachindex(t))
        setindex!(t, Threads.@spawn(setindex!(v, Model(), j)), k)
    end
    foreach(wait, t)
end
printinfo() = (th = map(Threads.nthreads, (:default, :interactive)); println("Settings> Threads=$th"))

end

# ✅ Status: 3==INFEASIBLE, 5==UNBOUNDED, 4==(3||5)
