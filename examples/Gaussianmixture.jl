f(x,μ,σ) = 1/sqrt(2 *π*σ^2) * exp(-(x - μ)^2 / (2σ^2))
μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
χ(x) = sum(w[i]*f(x,μ[i],σ[i]) for i in 1:length(w))

using PolyChaos
meas = Measure("myMeas",χ,(-∞,∞),false,Dict(:μ=>μ,:σ=>σ,:w=>w))

