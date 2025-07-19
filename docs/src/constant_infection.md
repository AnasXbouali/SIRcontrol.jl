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
function SIRocp(tf, S0, I0, Q, β, γ, umax)
	    ocp = @def begin
	        t ∈ [0, tf], time
	        x = (S, I, C) ∈ R^3, state
	        u ∈ R, control
	
	        S(0) == S0
	        I(0) == I0
	        C(0) == Q
	
	        C(tf) ≥ 0.0
	
	        0 ≤ u(t) ≤ umax
	
	        0.0 ≤ I(t) ≤ 1.0
	        0.0 ≤ S(t) ≤ 1.0
	
	        ẋ(t) == [-(1 - u(t)) * β * S(t) * I(t),
	                  (1 - u(t)) * β * S(t) * I(t) - γ * I(t),
	                -u(t)]
	
	        -S(tf) → min
	    end
	
	    sol = solve(ocp, :direct, :adnlp, :ipopt; 
	                    disc_method = :gauss_legendre_3, 
	                    grid_size=1000, 
                            tol=1e-9, 
	                    display=true)
            return sol
end
```

```@example main
tf   = 300
γ    = 0.2
S0   = 0.999
I0   = 0.001
Q    = 28
β    = [0.6,0.9]
umax = [0.9,  0.5]

sol1 = SIRocp(tf, S0, I0, Q, β[1], γ, umax[1]) 
sol2 = SIRocp(tf, S0, I0, Q, β[2], γ, umax[2])
```























## Case 1: $\beta = 0.6$ and $\bar{u} = 0.9$

```@example main
gr()

plot(
    t -> state(sol1)(t)[1], 0, tf,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol1)(t)[2], 0, tf,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf], [γ/β[1], γ/β[1]],
    label = L"$S_h$", lw = 2, ls = :dash, color = :black
)
plot!(
    t -> control(sol1)(t)[1], 0, tf,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
```

```@example main
t = 0:1e-4:50
t_sol1 = t[findfirst(i -> abs(control(sol1)(t[i])[1]) > 0.9*umax[1], 2:length(t))] # intervention time
nothing
```

## Case 1: Plot of $\log(S(t_c)) - \log(S(t_c+\Delta))$

```@example main
Δ1 = Q / umax[1]

tspan = (0.0, tf)

function control_(t, tc)
    if t ≥ tc && t ≤ tc + Δ1
        return umax[1]
    else
        return 0.0
    end
end

function sir!(du, u, p, t)
    S, I = u
    tc = p
    u_val = control_(t, tc)
    du[1] = -(1 - u_val) * β[1] * S * I
    du[2] = (1 - u_val) * β[1] * S * I - γ * I
end

tcs = 0.0:0.05:(tf-Δ1)
result_values = Float64[] 

for tc in tcs
   
    prob = ODEProblem(sir!, [S0, I0], tspan, (tc))
    sol = solve(prob, Tsit5(), abstol=1e-12, reltol=1e-12)
    if tc >= first(sol.t) && tc <= last(sol.t) && (tc + Δ1) >= first(sol.t) && (tc + Δ1) <= last(sol.t)
        S_tc = sol(tc)[1] 
        S_tcD = sol(tc + Δ1)[1] 
        
        val = log(S_tc) - log(S_tcD)
        push!(result_values, val) 
    else
        push!(result_values, NaN) 
    end
end
```

```@example main
plot(tcs, result_values;
    label = L"t_c \mapsto \log S(t_c) - \log S(t_c + \Delta)",
    lw = 2,
    grid = false)
plot!([t_sol1, t_sol1], [0, 0.1], linestyle = :dash, lw = 2, label = L"t_c^*")
```









## Case 2: $\beta = 0.9$ and $\bar{u} = 0.5$


```@example main
gr()

plot(
    t -> state(sol2)(t)[1], 0, tf,
    label = L"$S$", color = :darkblue, lw = 2,
    grid = false
)

plot!(
    t -> state(sol2)(t)[2], 0, tf,
    label = L"$I$", color = :darkgreen, lw = 2
)

plot!(
    [0, tf], [γ/β[2], γ/β[2]],
    label = L"$S_h$", lw = 2, ls = :dash, color = :black
)
plot!(
    t -> control(sol2)(t)[1], 0, tf,
    label = L"$u$", color = :darkred, lw = 2,
    grid = false
)
```

```@example main
t = 0:1e-4:tf
t_sol2 = t[findfirst(i -> abs(control(sol2)(t[i])[1]) >  umax[2], 2:length(t))] # intervention time
nothing
```

## Case 2: Plot of $\log(S(t_c)) - \log(S(t_c+\Delta))$


```@example main
Δ2 = Q / umax[2]

tspan = (0.0, tf)

function control__(t, tc)
    if t ≥ tc && t ≤ tc + Δ1
        return umax[2]
    else
        return 0.0
    end
end

function sir!(du, u, p, t)
    S, I = u
    tc = p
    u_val = control__(t, tc)
    du[1] = -(1 - u_val) * β[2] * S * I
    du[2] = (1 - u_val) * β[2] * S * I - γ * I
end

tcs = 0.0:0.05:(tf-Δ1)
result_values_ = Float64[] 

for tc in tcs
   
    prob_ = ODEProblem(sir!, [S0, I0], tspan, (tc))
    sol_ = solve(prob_, Tsit5(), abstol=1e-12, reltol=1e-12)
    if tc >= first(sol_.t) && tc <= last(sol_.t) && (tc + Δ1) >= first(sol_.t) && (tc + Δ1) <= last(sol_.t)
        S_tc = sol_(tc)[1] 
        S_tcD = sol_(tc + Δ1)[1] 
        
        val = log(S_tc) - log(S_tcD)
        push!(result_values_, val) 
    else
        push!(result_values_, NaN) 
    end
end
```

```@example main
plot(tcs, result_values_;
    label = L"t_c \mapsto \log S(t_c) - \log S(t_c + \Delta)",
    lw = 2,
    grid = false)
plot!([t_sol2, t_sol2], [0, 1.81], linestyle = :dash, lw = 2, label = L"t_c^*")
```


