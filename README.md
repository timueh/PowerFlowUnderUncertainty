A case study that compares probabilistic power flow, chance-constrained AC optimal power flow, and chance-constrained DC optimal power flow.

More information to follow soon.

The code implements methods from the following papers:

- [Chance-constrained AC-OPF](https://ieeexplore.ieee.org/document/8719988)
- [Chance-constrained DC-OPF](https://www.sciencedirect.com/science/article/pii/S235246771830105X) (open access)

Polynomial chaos is used to propagate uncertainties, for which we use [`PolyChaos.jl`](https://github.com/timueh/PolyChaos.jl).
