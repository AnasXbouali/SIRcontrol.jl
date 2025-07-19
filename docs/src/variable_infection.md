```@example main
using OptimalControl
using NLPModelsIpopt
using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Roots
```


```@example main
function beta(t, τ, β1, β2, k=300)
    return β1 + (β2 - β1) / (1 + exp(-k * (t - τ))) 
end
```

```@example main
tf_   = 300
τ     = 50
γ_    = 0.15
β11   = 0.4
β12   = 0.2
S0_   = 0.999
I0_   = 0.001
Q_    = 35
β21   = 0.7
β22   = 0.4
umax_ = 0.8
k = 300         # regularization parameter

plot(t -> beta(t, τ, β11, β12, k), 0, tf_, label = L"$\beta$", color = :black, lw = 2, grid = false)
```

```@example main
function SIRocp2beta(T, S0, I0, Q, τ, β1, β2, γ, umax,k)
	
	ocp = @def begin
	        t ∈ [0, T], time
	        x = (S, I, C) ∈ R^3, state
	        u ∈ R, control
	
	        S(0) == S0
	        I(0) == I0
	        C(0) == Q
	
	        C(T) ≥ 0.0
	
	        0 ≤ u(t) ≤ umax
	
	        0.0 ≤ I(t) ≤ 1.0
	        0.0 ≤ S(t) ≤ 1.0
	
	        ẋ(t) == [-(1 - u(t)) * beta(t, τ, β1, β2, k) * S(t) * I(t),
	                  (1 - u(t)) * beta(t, τ, β1, β2, k) * S(t) * I(t) - γ * I(t),
	                -u(t)]
	
	        -S(T) → min
	end
	
	sol = nothing
	for N in [50, 100, 200,300, 400, 500,600,700,800,900,1000]
	    sol = solve(ocp, :direct, :adnlp, :ipopt; 
	                    disc_method = :gauss_legendre_3, 
	                    grid_size=N, 
	                    init=sol,
                        tol=1e-8, 
	                    display=true)
	end
	return sol
end
```

```@example main
sol1_ = SIRocp2beta(tf_, S0_, I0_, Q_, τ, β11, β12, γ_, umax_, k)
sol2_ = SIRocp2beta(tf_, S0_, I0_, Q_, τ, β21, β22, γ_, umax_, k)
```

## case 1: $\beta_1 > \beta_2$ (two interventions)


```@example main
gr()

plot(
    t -> state(sol1_)(t)[1], 0, tf_,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol1_)(t)[2], 0, tf_,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf_], [γ_/ β11, γ_/ β11],
    label = L"$S_h^1$", lw = 2, ls = :dash, color = :black
)
plot!(
    [0, tf_], [γ_/ β12, γ_/ β12],
    label = L"$S_h^2$", lw = 2, ls = :dash, color = :purple
)
plot!(
    t -> control(sol1_)(t)[1], 0, tf_,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
vline!([τ], lw=2, color = "black", linestyle = :dash, label = L"$\tau$")
```

## case 1bis: $\beta_1 > \beta_2$ (one interventions)

```@example main
gr()

plot(
    t -> state(sol2_)(t)[1], 0, tf_,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol2_)(t)[2], 0, tf_,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf_], [γ_/ β21, γ_/ β21],
    label = L"$S_h^1$", lw = 2, ls = :dash, color = :black
)
plot!(
    [0, tf_], [γ_/ β22, γ_/ β22],
    label = L"$S_h^2$", lw = 2, ls = :dash, color = :purple
)
plot!(
    t -> control(sol2_)(t)[1], 0, tf_,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
vline!([τ], lw=2, color = "black", linestyle = :dash, label = L"$\tau$")
```

```@example main
sol3_ = SIRocp2beta(tf_, S0_, I0_, Q_, 150, 0.4, 0.8, γ_, umax_, k)
sol4_ = SIRocp2beta(tf_, S0_, I0_, Q_, 100, 0.2, 0.8, γ_, umax_, k)
```

## case 2: $\beta_1 < \beta_2$ (two interventions)


```@example main
gr()

plot(
    t -> state(sol3_)(t)[1], 0, tf_,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol3_)(t)[2], 0, tf_,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf_], [γ_/ 0.4, γ_/ 0.4],
    label = L"$S_h^1$", lw = 2, ls = :dash, color = :black
)
plot!(
    [0, tf_], [γ_/ 0.8, γ_/ 0.8],
    label = L"$S_h^2$", lw = 2, ls = :dash, color = :purple
)
plot!(
    t -> control(sol3_)(t)[1], 0, tf_,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
vline!([150], lw=2, color = "black", linestyle = :dash, label = L"$\tau$")

```

## case 2bis: $\beta_1 < \beta_2$ (one interventions)

```@example main
gr()

plot(
    t -> state(sol4_)(t)[1], 0, tf_,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol4_)(t)[2], 0, tf_,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf_], [γ_/ 0.2, γ_/ 0.2],
    label = L"$S_h^1$", lw = 2, ls = :dash, color = :black
)
plot!(
    [0, tf_], [γ_/ 0.8, γ_/ 0.8],
    label = L"$S_h^2$", lw = 2, ls = :dash, color = :purple
)
plot!(
    t -> control(sol4_)(t)[1], 0, tf_,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
vline!([100], lw=2, color = "black", linestyle = :dash, label = L"$\tau$")
```