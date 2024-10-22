import Symbolics
using GalaxyProfiles

# Define variables
Symbolics.@variables ρ0, rs, r, rc, n
dr = Symbolics.Differential(r)
# Symbolic ∇ρ
Symbolics.expand_derivatives(dr(ρ(CoreNFW(NFW(ρ0, rs), rc, n), r))) |> println
# Convert floats to rationals with Rational(x) for consistent promotion;
# Mathematica gives better simplification but more work to copy over

