export  sampleFromGaussianMixture,
        setupUncertainty,
        generateSamples

function ρ_gauss(x,μ,σ)
    1/sqrt(2*π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
end

function setupUncertainty(μ::Vector{},σ::Vector{},w::Vector{},n::Int,deg::Int)
    @assert length(μ) == length(σ) == length(w) "inconsistent lengths of μ and σ"
    ρ(x) = sum( w[i]*ρ_gauss(x,μ[i],σ[i]) for i in 1:length(w) )
    meas = Measure("my_GaussMixture", ρ, (-Inf,Inf), false, Dict(:μ=>μ,:σ=>σ,:w=>w)) # build measure
    op = OrthoPoly("my_op",deg,meas;Nquad = 100,Nrec = 5*deg, discretization=stieltjes) # construct orthogonal polynomial
    opq = OrthoPolyQ(op)
    showbasis(op,digits=2) # in case you wondered

    pd = zeros(n,deg+1)
    pd[1, [1,2]] = calculateAffinePCE(opq)
    pd[2, 1] = 1.2
    qd = 0.85 * copy(pd)

    return Dict(:opq=>opq,
                :T2=>Tensor(2,opq),
                :T3=>Tensor(3,opq),
                :pd=>pd,
                :qd=>qd,
                :dim=>size(pd,2))
end

function sampleFromGaussianMixture(n::Int,μ::Vector{},σ::Vector{},w::Vector{})
    X = Float64[]
    for i in 1:n
        k = findfirst(x -> x > rand(), cumsum(w))
        push!(X, μ[k] + σ[k]*randn())
    end
    return X
end

function generateSamples(x,d_in::Dict,unc::Dict)
    Φ = evaluate(x,unc[:opq])
    d_out = Dict{Symbol,Matrix{Float64}}()
    for (key, value) in d_in
        d_out[key] = (Φ*value')'
    end
    d_out[:i] = @. sqrt( d_out[:i_re]^2 + d_out[:i_im]^2 )
    d_out[:v] = @. sqrt( d_out[:e]^2 + d_out[:f]^2)
    return d_out
end
