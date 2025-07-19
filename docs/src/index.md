# 🧪 On the problem of minimizing the epidemic final size for SIR model via social distancing

## 📘 Introduction

### **Context & Motivation**

We revisit the problem of minimizing the epidemic final size in the **SIR model** through **social distancing interventions**.

Traditionally, this problem assumes a fixed **interval structure** for the timing of interventions. In contrast, we investigate a more flexible approach: **$L^1$-constrained controls** that limit total intervention effort, rather than when and how it's applied.

Key insight:

> Even with this generalization, the optimal control still occurs over a **single time interval**—indicating no gain from fragmenting interventions across multiple periods.

However, if the **infection rate \$\beta(t)\$ changes over time**, then **splitting interventions** into **two disjoint periods** may become optimal.

---

### **Our Contributions**

#### 🔍 Theoretical Results

1. **Constant Transmission Rate**

   > Optimal control is supported on **one continuous interval**.
   > Concentrated efforts outperform spread-out strategies.

2. **Piecewise Constant Transmission Rate**

   > When \$\beta(t)\$ changes once (e.g., from \$\beta\_1\$ to \$\beta\_2\$), the optimal policy may involve **at most two separate interventions**.

#### 💻 Numerical Simulations

We support our theoretical claims with simulations, covering:

* Constant \$\beta\$
* Single-change \$\beta\$ (piecewise constant)

These validate the analytical results and quantify the differences between intervention strategies.

---

## ⚙️ The SIR Model with Control 

### **System Dynamics**

We use the classical **SIR model** with a control variable \$u(t)\$ (representing social distancing):

```math
\begin{cases}
\dot{S}(t) = -(1-u(t))\beta(t) S(t)I(t) \\
\dot{I}(t) = (1-u(t))\beta(t) S(t)I(t) - \gamma I(t)
\end{cases}
```

* \$S(t)\$, \$I(t)\$: Susceptible and infected population fractions
* \$u(t) \in \[0, \bar{u}]\$: Intervention intensity
* \$\beta(t)\$: Transmission rate
* \$\gamma\$: Recovery rate

> The control \$u(t)\$ modifies the infection term, reducing contacts based on intervention intensity.

---

### **Key Parameters**

* **Transmission Rate**: \$\beta(t) > 0\$ (time-varying)
* **Recovery Rate**: \$\gamma > 0\$
* **Reproduction Number**: \$R\_0(t) = \beta(t)/\gamma\$

---

### **Transmission Scenarios**

We analyze two transmission models:

1. **Constant \$\beta\$**

```math
   \beta(t) = \beta_0 \quad \forall t \geq 0
```

2. **Piecewise Constant \$\beta\$**

```math
   \beta(t) = \begin{cases}
   \beta_1 & \text{if } 0 \leq t < T_c \\
   \beta_2 & \text{if } t \geq T_c
   \end{cases}
```

   where \$T\_c\$ is the time of the transmission change.

---

## 🎯 Optimal Control Formulation

### **Objective**

Our goal: **Minimize the final epidemic size** = maximize \$S(\infty)\$.

```math
J(u) = S(\infty) = \lim_{t \to \infty} S(t)
```

This represents preserving the **largest possible fraction** of the susceptible population.

---

### **Constraints**

We impose an **\$L^1\$ budget constraint** on the intervention:

```math
\|u(\cdot)\|_{L^1} = \int_0^{\infty} u(t) \, dt \leq K
```

* \$K > 0\$: Total control budget
* \$\bar{u}\$: Maximum intensity at any time

This approach generalizes older models that assume rigid time intervals for control.

---

### **Admissible Control Set**

We define the space of allowed controls:

```math
\mathcal{U} = \{u : [0,\infty) \to [0,\bar{u}] \,|\, u \text{ is measurable and } \|u\|_{L^1} \leq K \}
```
