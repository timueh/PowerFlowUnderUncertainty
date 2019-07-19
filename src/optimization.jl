export  addVariables!,
        addConsistencyConstraints!,
        addPowerFlow!,
        addSlackBus!,
        addPVBus!,
        addCost!,
        addInitialCondition!,
        addCore!

function addVariables!(opf::Model,sys::Dict,unc::Dict)
    N, Ng = sys[:N], sys[:Ng]
    K = size(unc[:pd],2)
    @variable(opf, pg[1:Ng,1:K])
    @variable(opf, qg[1:Ng,1:K])
    @variable(opf, e[1:N,1:K])
    @variable(opf, f[1:N,1:K])
end

function addConsistencyConstraints!(opf::Model)
    e, f = opf[:e], opf[:f]
    @constraint(opf, 0.7 .<= e[:,1] .<= 1.3)
    @constraint(opf, -0.3 .<= f[:,1] .<= 0.3)
end

function addPowerFlow!(opf::Model,sys::Dict,unc::Dict)
    N, Ng = sys[:N], sys[:Ng]
    K = size(unc[:pd],2)
    e, f, pg, qg = opf[:e], opf[:f], opf[:pg], opf[:qg]
    G, B = real(sys[:Ybus]), imag(sys[:Ybus])
    p, q = sys[:Cp]*pg - sys[:Cd]*unc[:pd], sys[:Cp]*qg - sys[:Cd]*unc[:qd]
    @constraint(opf, pfP[i=1:N,k=1:K], p[i,k] == sum( (G[i,j]*(e[i,l1]*e[j,l2]+f[i,l1]*f[j,l2]) + B[i,j]*(f[i,l1]*e[j,l2]-e[i,l1]*f[j,l2]) )*unc[:T3].get([l1-1,l2-1,k-1])/unc[:T2].get([k-1,k-1]) for l2 in 1:K, l1 in 1:K, j in 1:N) )
    @constraint(opf, pfQ[i=1:N,k=1:K], q[i,k] == sum( (G[i,j]*(f[i,l1]*e[j,l2]-e[i,l1]*f[j,l2]) - B[i,j]*(e[i,l1]*e[j,l2]+f[i,l1]*f[j,l2]) )*unc[:T3].get([l1-1,l2-1,k-1])/unc[:T2].get([k-1,k-1]) for l2 in 1:K, l1 in 1:K, j in 1:N) )
end

function addSlackBus!(opf::Model)
    e, f = opf[:e], opf[:f]
    K = size(e,2)
    @constraint(opf, slack_e, e[1,:] .== [1; zeros(K-1)])
    @constraint(opf, slack_f, f[1,:] .== zeros(K))
end

function addPVBus!(opf::Model,sys::Dict,unc::Dict)
    e, f, pg = opf[:e], opf[:f], opf[:pg]
    K, bus, gen = size(e,2), 4, 2
    v = [1.04^2; zeros(K-1)]
    @constraint(opf, PVBus_V[k=1:K], sum((e[bus,l1]*e[bus,l2]+f[bus,l1]*f[bus,l2])*unc[:T3].get([l1-1,l2-1,k-1])/unc[:T2].get([k-1,k-1]) for l1=1:K, l2=1:K) == v[k] )
    @constraint(opf, PVBus_P, pg[gen,:] .== [ 0.84; zeros(K-1) ])
end

function addCost!(opf::Model,sys::Dict,unc::Dict)
    p = opf[:pg]
    Ng, K = size(p)
    @objective(opf, Min, sum( p[i,1] * sys[:costlin][i] + sys[:costquad][i]*sum(p[i,l]^2 * unc[:T2].get([l-1,l-1]) for l in 1:K) for i in 1:Ng) )
end

function addInitialCondition!(opf::Model,sys::Dict,unc::Dict)
    e, f, pg, qg = opf[:e], opf[:f], opf[:pg], opf[:qg]
    N, K = size(e)
    # [ set_start_value(pg[i], -unc[:pd][i]) for i in 1:length(unc[:pd]) ]
    # [ set_start_value(qg[i], -unc[:qd][i]) for i in 1:length(unc[:qd]) ]
    [ set_start_value(pg_, 0) for pg_ in pg ]
    [ set_start_value(qg_, 0) for qg_ in qg ]
    [ set_start_value(e_,1) for e_ in e[:,1] ]
    [ set_start_value(e_,0) for e_ in e[:,2:end] ]
    [ set_start_value(f_,0) for f_ in f ]
end

function addCore!(mod::Model,sys::Dict,unc::Dict)
    addVariables!(mod,sys,unc)
    addInitialCondition!(mod,sys,unc)
    addPowerFlow!(mod,sys,unc)
    addSlackBus!(mod)
end
