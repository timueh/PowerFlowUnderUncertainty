export  addVariables!,
        addConsistencyConstraints!,
        addPowerFlow!,
        addSlackBus!,
        addPVBus!,
        addCost!,
        addInitialCondition!,
        addCore!,
        addVariablesDC!,
        addPowerFlowDC!,
        addCostDC!,
        addCoreDC!,
        addConstraintsGenDC!,
        addConstraintsLineFlowDC!,
        buildCostSOC

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

#########################################################################
## DC Equations ##

function addVariablesDC!(mod::Model,sys::Dict,unc::Dict)
    N, Ng = sys[:N], sys[:Ng]
    K = size(unc[:pd],2)
    @variable(mod, pg[i in 1:Ng, j in 1:K], base_name="p")
end

function addPowerFlowDC!(opf::Model,sys::Dict,unc::Dict)
    Nd, Ng = sys[:Nd], sys[:Ng]
    p, d = opf[:pg], unc[:pd]
    K = size(d,2)
    @assert size(d,1) == Nd "Detected inconsistency in number of demands."
    @constraint(opf, energy_balance[j in 1:K], sum(p[i,j] for i in 1:Ng) - sum(d[i,j] for i in 1:Nd) == 0)
end

function addConstraintsGenDC!(mod::Model,sys::Dict,unc::Dict)
    lb, λlb = sys[:con][:pg][:lb], sys[:con][:pg][:λ][:lb]
    ub, λub = sys[:con][:pg][:ub], sys[:con][:pg][:λ][:ub]
    p = mod[:pg]
    Ng, K = size(p)
    mop = typeof(unc[:opq]) == OrthoPolyQ ? MultiOrthoPoly([unc[:opq]],unc[:opq].op.deg) : unc[:opq]
    @constraint(mod, con_pmax[i in 1:Ng], [1/λub[i]*(ub[i] - mean(p[i,:],mop)); buildSOC(p[i,:],mop)] in SecondOrderCone())
    @constraint(mod, con_pmin[i in 1:Ng], [1/λlb[i]*(mean(p[i,:],mop) - lb[i]); buildSOC(p[i,:],mop)] in SecondOrderCone())
end

function addConstraintsLineFlowDC!(mod::Model,sys::Dict,unc::Dict)
    lb, λlb = sys[:con][:pl][:lb], sys[:con][:pl][:λ][:lb]
    ub, λub = sys[:con][:pl][:ub], sys[:con][:pl][:λ][:ub]
    p, d = mod[:pg], unc[:pd]
    Ng, K = size(p)
    pl = sys[:ptdf]*(sys[:Cp]*p + sys[:Cd]*d)
    Nline = size(pl,1)
    mop = typeof(unc[:opq]) == OrthoPolyQ ? MultiOrthoPoly([unc[:opq]],unc[:opq].op.deg) : unc[:opq]
    @constraint(mod, con_plmax[i in 1:Nline], [1/λub[i]*(ub[i] - mean(pl[i,:],mop)); buildSOC(pl[i,:],mop)] in SecondOrderCone())
    @constraint(mod, con_plmin[i in 1:Nline], [1/λlb[i]*(mean(pl[i,:],mop) - lb[i]); buildSOC(pl[i,:],mop)] in SecondOrderCone())
end

function addCoreDC!(mod::Model,sys::Dict,unc::Dict)
    addVariablesDC!(mod,sys,unc)
    addPowerFlowDC!(mod,sys,unc)
end

function buildSOC(x::Vector,mop::MultiOrthoPoly)
    t = [ sqrt(Tensor(2,mop).get([i,i])) for i in 0:mop.dim-1 ]
    (t.*x)[2:end]
end

function buildCostSOC_perbus(bus::Int,lin::Real,quad::Real,T2::Tensor,K::Int)
    [ 0.5*lin; zeros(K-1) ], quad*[ T2.get([k-1,k-1]) for k in 1:K ]
end

function buildCostSOC(lin::Vector,quad::Vector,T2::Tensor,K::Int)
    @assert length(quad) == length(lin) "inconsistent length of cost coefficients"
    res = [ buildCostSOC_perbus(i,lin[i],quad[i],T2,K) for i in 1:length(lin) ]
    vcat([r[1] for r in res]...), vcat([r[2] for r in res]...) # unpack res
end

function addCostDC!(opf::Model,sys::Dict,unc::Dict)
    p = opf[:pg]
    pvec = vec(p')
    Ng, K = size(p)
    lin, quad = buildCostSOC(sys[:costlin],sys[:costquad],Tensor(2,unc[:opq]),K)
    @variable(opf, t)
    @constraint(opf, [t; sqrt.(quad) .* pvec + (1 ./ sqrt.(quad)).*lin ] in SecondOrderCone())
    @objective(opf, Min, t)
end
