


```julia
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
