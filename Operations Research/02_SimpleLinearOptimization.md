# Simple Linear Optimization

````julia
# optimization framework
using JuMP;

# solvers
using CPLEX; using Gurobi; using Cbc;
````





$
max \; x_1 + 2x_2 + 5x_3
$

subject to,

$
\begin{aligned}
\ -x_1 + x_2 + 3x_3  \le -5 \\
\ x_1 3x_2 - 7x_3 \le 10 \\
\ 0 \le x_1 \le 10 \\
\ x_2 \ge 0 \\
\ x_3 \ge 0
\end{aligned}
$

````julia
# optimization model
model = Model(Gurobi.Optimizer)
````


````
Academic license - for non-commercial use only
````



````julia

# variables
@variable(model, 0 <= x1 <= 10)
@variable(model, x2 >= 0)
@variable(model, x3 >= 0)

# objective
@objective(model, Max, x1 + 2x2 + 5x3)

@constraint(model, constraint1, -x1 + x2 + 3x3 <= -5)
@constraint(model, constraint2, x1 + 3x2 - 7x3 <= 10)

# take a peek
model

optimize!(model)
````


````
Academic license - for non-commercial use only
Gurobi Optimizer version 9.0.2 build v9.0.2rc0 (win64)
Optimize a model with 2 rows, 3 columns and 6 nonzeros
Model fingerprint: 0x040b61ef
Coefficient statistics:
  Matrix range     [1e+00, 7e+00]
  Objective range  [1e+00, 5e+00]
  Bounds range     [1e+01, 1e+01]
  RHS range        [5e+00, 1e+01]
Presolve time: 0.00s
Presolved: 2 rows, 3 columns, 6 nonzeros

Iteration    Objective       Primal Inf.    Dual Inf.      Time
       0    9.0000000e+30   1.250000e+30   9.000000e+00      0s
       2    1.9062500e+01   0.000000e+00   0.000000e+00      0s

Solved in 2 iterations and 0.00 seconds
Optimal objective  1.906250000e+01
````



````julia

values = [ JuMP.value(x1), JuMP.value(x2), JuMP.value(x3) ]

display(values)
````


````
3-element Array{Float64,1}:
 10.0
  2.1875
  0.9375
````



````julia

JuMP.dual(constraint1)
JuMP.dual(constraint2)
````


````
-0.06249999999999989
````


