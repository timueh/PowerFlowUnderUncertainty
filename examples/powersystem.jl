function setupPowerSystem()
    # incidence matrix
    A = [ -1 1 0 0; -1 0 1 0; -1 0 0 1 ; 0 1 -1 0; 0 0 -1 1]
    Nline, N = size(A,1), size(A,2)

    # vector of branch parameters
    Rbr, Xbr = [0.01008
    0.00744
    0.00744
    0.01272
    0.01004], [0.0504
    0.0372
    0.0372
    0.0636
    0.0601]
    @assert length(Rbr) == length(Xbr) == Nline "inconsistent branch parameters"
    Ybr = diagm(0 => 1 ./ (Rbr + im*Xbr))
    Bsh = zeros(size(Ybr))
    # bus admittance matrix
    Ybus = A'*Ybr*A

    # PTDF matrix
    Bbr = -diagm(0 => 1 ./ Xbr)
    Bbus = A'*Bbr*A
    Ψ = [ zeros(Nline)  -Bbr*A[:,2:end]*inv(A[:,2:end]'*Bbr*A[:,2:end]) ]

    # book-keeping
    Cp, Cd = [1 0; 0 0; 0 0; 0 1], [0 0; 1 0; 0 1; 0 0 ]
    Ng, Nd = size(Cp,2), size(Cd,2)

    # cost
    costquad, costlin = [1000, 5000], ones(Ng)
    # λp, λl = 1.6*ones(Ng), 1.6*ones(Nline) # lambdas for chance constraint reformulations
    # pmax, pmin = 10*ones(Ng), zeros(Ng) # engineering limits
    # plmax, plmin = 10*ones(Nline), -10*ones(Nline) # engineering limits
    con = Dict(:pg => Dict(
                    :λ => Dict(
                        :lb => 1.6*ones(Ng),
                        :ub => 1.6*ones(Ng)),
                    :lb => zeros(Ng),
                    :ub => [4.0, 10.0] ),
               :qg => Dict(
                    :λ => Dict(
                        :lb => 1.6*ones(Ng),
                        :ub => 2.6*ones(Ng)),
                    :lb => zeros(Ng),
                    :ub => [10.0, 1.2] ),
               :pl => Dict(
                    :λ => Dict(
                        :lb => 1.6*ones(Nline),
                        :ub => 2.0*ones(Nline)),
                    :lb => -10*ones(Nline),
                    :ub =>  [3.5, 10, 10, 2.0, 1.30] ),
                :v => Dict(
                     :λ => Dict(
                         :lb => 2.5*ones(N),
                         :ub => 1.6*ones(N)),
                     :lb => [0.9, 0.81, 0.95, 0.8],
                     :ub => 1.005*ones(N) )
                        )

    return Dict(:Ybus=>Ybus,
                :Bbus=>Bbus,
                :Ybr=>Ybr,
                :Bsh=>Bsh,
                :A=>A,
                :ptdf=>Ψ,
                :Cp=>Cp,
                :Cd=>Cd,
                :Ng=>Ng,
                :Nd=>Nd,
                :N=>N,
                :Nline=>Nline,
                :costquad=>costquad,
                :costlin=>costlin,
                :con=>con
                )
end
##
