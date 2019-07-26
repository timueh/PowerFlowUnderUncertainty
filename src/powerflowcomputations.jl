export  getGridState,
        computeLineFlows,
        computeLineCurrents,
        getGridStateDC

function getGridState(mod::Model,sys::Dict,unc::Dict)
    d = Dict{Symbol,Matrix{Float64}}()
    for (key, val) in mod.obj_dict
        typeof(val) == Matrix{VariableRef} ? (d[key] = @. value(val)) : nothing
    end
    d[:pd], d[:qd] = unc[:pd], unc[:qd]
    merge!(d,computeLineFlows(mod,sys,unc))
    merge!(d,computeLineCurrents(mod,sys,unc))
end

function getGridStateDC(mod::Model,sys::Dict,unc::Dict)
    d = Dict{Symbol,Matrix{Float64}}()
    for (key, val) in mod.obj_dict
        typeof(val) == Matrix{VariableRef} ? (d[key] = @. value(val)) : nothing
    end
    merge!(d,computeLineFlowsDC(mod,sys,unc))
    merge!(d,computeAnglesDC(mod,sys,unc))
end

function computeLineCurrents(mod::Model,grid::Dict,unc::Dict)
    A, ybr = grid[:A], grid[:Ybr]
    gbr, bbr, bsh = real(ybr), imag(ybr), grid[:Bsh]
    M, N, L = grid[:Nline], grid[:N], unc[:dim]
    E, F = value.(mod[:e]), value.(mod[:f])
    # line power flow computations
    fil_real, fil_imag = zeros(M,L), zeros(M,L)
    for l in 1:M
        i, j = i, j = getNodesForLine(A,l)
        for k in 1:L
            fil_real[l,k] = gbr[l,l]*(E[i,k]-E[j,k]) - bbr[l,l]*(F[i,k]-F[j,k])
            fil_imag[l,k] = bbr[l,l]*(E[i,k]-E[j,k]) + gbr[l,l]*(F[i,k]-F[j,k])
        end
    end
    return Dict(:i_re=>fil_real,
                :i_im=>fil_imag)
end

function computeLineFlows(mod::Model,grid::Dict,unc::Dict)
    A, ybr, T3, T2 = grid[:A], grid[:Ybr], unc[:T3], unc[:T2]
    gbr, bbr, bsh = real(ybr), imag(ybr), grid[:Bsh]
    M, N, L = grid[:Nline], grid[:N], unc[:dim]
    P, Q, E, F = value.(mod[:pg]), value.(mod[:qg]), value.(mod[:e]), value.(mod[:f])
    # line power flow computations
    fpl_t, fpl_f, fql_t, fql_f = zeros(M,L), zeros(M,L), zeros(M,L), zeros(M,L)
    for l in 1:M
        i, j = getNodesForLine(A,l)
        for k in 1:L
            fpl_t[l,k] = sum( sum( ( gbr[l,l]*(E[i,l1]*E[i,l2]-E[i,l1]*E[j,l2]+F[i,l1]*F[i,l2]-F[i,l1]*F[j,l2]) + bbr[l,l]*(E[i,l1]*F[j,l2]-F[i,l1]*E[j,l2]) )*T3.get([l1-1,l2-1,k-1])/T2.get([k-1,k-1]) for l2=1:L) for l1=1:L)
            fpl_f[l,k] = sum( sum( ( gbr[l,l]*(E[j,l1]*E[j,l2]-E[i,l1]*E[j,l2]+F[j,l1]*F[j,l2]-F[i,l1]*F[j,l2]) - bbr[l,l]*(E[i,l1]*F[j,l2]-F[i,l1]*E[j,l2]) )*T3.get([l1-1,l2-1,k-1])/T2.get([k-1,k-1]) for l2=1:L) for l1=1:L)
            fql_t[l,k] = sum( sum( ( gbr[l,l]*(E[i,l1]*F[j,l2]-F[i,l1]*E[j,l2]) - (bbr[l,l]+bsh[l,l])*(E[i,l1]*E[i,l2] + F[i,l1]*F[i,l2]) + bbr[l,l]*(E[i,l1]*E[j,l2]+F[i,l1]*F[j,l2]) )*T3.get([l1-1,l2-1,k-1])/T2.get([k-1,k-1]) for l2=1:L) for l1=1:L)
            fql_f[l,k] = sum( sum( (-gbr[l,l]*(E[i,l1]*F[j,l2]-F[i,l1]*E[j,l2]) - (bbr[l,l]+bsh[l,l])*(E[j,l1]*E[j,l2] + F[j,l1]*F[j,l2]) + bbr[l,l]*(E[i,l1]*E[j,l2]+F[i,l1]*F[j,l2]) )*T3.get([l1-1,l2-1,k-1])/T2.get([k-1,k-1]) for l2=1:L) for l1=1:L)
        end
    end
    return Dict(:pl_t=>fpl_t,
                :ql_t=>fql_t,
                :pl_f=>fpl_f,
                :ql_f=>fql_f )
end

function computeLineFlowsDC(mod::Model,grid::Dict,unc::Dict)
    p, d = value.(mod[:pg]), unc[:pd]
    Dict(:pl_t => grid[:ptdf] * ( grid[:Cp]*p - grid[:Cd]*d ))
end

function computeAnglesDC(mod::Model,grid::Dict,unc::Dict)
    p, d, Bbus = value.(mod[:pg]), unc[:pd], grid[:Bbus]
    Bbusinv, pnet = inv(Bbus[2:end,2:end]), grid[:Cp]*p - grid[:Cd]*d
    Dict(:θ => [ zeros(1,size(p,2)); -Bbusinv * pnet[2:end,:] ])
end

# given the incidence matrix of a graph, find the nodes i and j that line l makes up
# -1 ==> line begins
#  1 ==> line ends
function getNodesForLine(A::Matrix,l::Int)
    findfirst(isequal(-1),A[l,:]), findfirst(isequal(1),A[l,:])
end
