---
title: Introduction
author: 'Brandon Moretz'
date: '15th May 2020'
---  
  
  
#  Introduction
  
  
```julia
using JuMP; using Clp;
  
using GLPKMathProgInterface
```
  
Example Solver
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?min%20&#x5C;sum_{(i,j)}%20c_{ij}x{ij}"/></p>  
  
  
subject to
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?&#x5C;sum_{i,j%20&#x5C;in%20A}%20x_{ij}%20-%20&#x5C;sum_{(j,i)&#x5C;in%20A}%20x_{ij}%20=%20b_i%20&#x5C;;%20&#x5C;forall_i%20&#x5C;in%20N"/></p>  
  
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?0%20&#x5C;le%20x_{ij}%20&#x5C;le%201%20&#x5C;space%20&#x5C;forall(i,%20j)%20&#x5C;in%20A"/></p>  
  
  
  
```julia
model = Model(Clp.Optimizer)
  
@variable(model, 0 <= x <= 40)
@variable(model, y <= 0)
@variable(model, z <= 0)
  
@objective(model, Max, x + y + z)
  
@constraint(model, const1, -x + y + z <= 20)
@constraint(model, const2, x + 3y + z <= 30)
  
display(model)
  
optimize!(model)
  
results = [JuMP.value(x), JuMP.value(y), JuMP.value(z) ]
  
display(results)
```
  