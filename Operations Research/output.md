---
title: Introduction
author: 'Brandon Moretz'
date: '15th May 2020'
---  
  
  
#  Simple Linear Optimization
  
  
```julia
# optimization framework
using JuMP;
  
# solvers
using CPLEX; using Gurobi; using Cbc;
```
  
##  Linear Programs
  
  
A straight forward approach to linear programming problems:
  
Given,
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?max%20&#x5C;;%20x_1%20+%202x_2%20+%205x_3"/></p>  
  
  
subject to,
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;begin{aligned}&#x5C;%20-x_1%20+%20x_2%20+%203x_3%20%20&#x5C;le%20-5%20&#x5C;&#x5C;&#x5C;%20x_1%203x_2%20-%207x_3%20&#x5C;le%2010%20&#x5C;&#x5C;&#x5C;%200%20&#x5C;le%20x_1%20&#x5C;le%2010%20&#x5C;&#x5C;&#x5C;%20x_2%20&#x5C;ge%200%20&#x5C;&#x5C;&#x5C;%20x_3%20&#x5C;ge%200&#x5C;end{aligned}"/>
  
In Julia:
  
```julia
# optimization model
model = Model(Gurobi.Optimizer)
  
# variables
@variable(model, 0 <= x1 <= 10)
@variable(model, x2 >= 0)
@variable(model, x3 >= 0)
  
# objective
@objective(model, Max, x1 + 2x2 + 5x3)
  
@constraint(model, constraint1, -x1 + x2 + 3x3 <= -5)
@constraint(model, constraint2, x1 + 3x2 - 7x3 <= 10)
  
# take a peek
display(model)
  
optimize!(model)
  
values = [ JuMP.value(x1), JuMP.value(x2), JuMP.value(x3) ]
  
display(values)
  
JuMP.dual(constraint1)
JuMP.dual(constraint2)
```
  
##  Alternative
  
  
```julia
# model
m2 = Model(Gurobi.Optimizer)
  
# non-zero constraints
@variable(m2, x[1:3] >= 0)
  
```
  
<img src="https://latex.codecogs.com/gif.latex?max%20&#x5C;;%20&#x5C;sum_{i=1}^3%20c_i%20x_i"/>
  
```julia
# coefficents
c = [1; 2; 5]
@objective(m2, Max, sum( c[i]*x[i] for i = 1:3))
```
  
Matrix Notation: <img src="https://latex.codecogs.com/gif.latex?Ax%20&#x5C;le%20b"/>
  
```julia
A = [-1 1 3;
     1  3 -7]
  
b = [-5; 10]
  
@constraint(m2, # model
  constraint[j=1:2], # num rows
  sum( A[j,i] * x[i] for i=1:3) <= b[j] )
  
# boundary constraint
@constraint(m2, bound, x[1] <= 10)
```
  
Solve it
  
```julia
optimize!(m2)
```
  
Results
  
```julia
values = [ JuMP.value(x1), JuMP.value(x2), JuMP.value(x3) ]
  
display(values)
  
JuMP.dual(constraint1)
JuMP.dual(constraint2)
```
  
##  Yet Another Way
  
  
```julia
  
m3 = Model(Gurobi.Optimizer)
  
c = [ 1; 2; 5]
A = [-1  1  3;
      1  3 -7]
b = [-5; 10]
  
index_x = 1:3
index_constraints = 1:2
  
@variable(m3, x[index_x] >= 0)
  
@objective(m3, Max, sum( c[i] * x[i] for i in index_x) )
  
@constraint(m3, constraint[j in index_constraints],
                sum( A[j, i] * x[i] for i in index_x) <= b[j] )
  
@constraint(m3, bound, x[1] <= 10)
  
optimize!(m3)
  
display(m3)
  
println("Optimal Solutions:")
for i in index_x
  println("x[$i] = ", value(x[i]))
end
  
println("Dual Variables:")
for j in index_constraints
  println("dual[$j] = ", dual(constraint[j]))
end
  
```
  
###  Mixed Integer Liner Programming
  
  
<img src="https://latex.codecogs.com/gif.latex?max%20&#x5C;;%20x_1%20+%202x_2%20+%205x_3"/>
  
subject to,
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;begin{aligned}&#x5C;%20-x_1%20+%20x_2%20+%203x_3%20%20&#x5C;le%20-5%20&#x5C;&#x5C;&#x5C;%20x_1%203x_2%20-%207x_3%20&#x5C;le%2010%20&#x5C;&#x5C;&#x5C;%200%20&#x5C;le%20x_1%20&#x5C;le%2010%20&#x5C;&#x5C;&#x5C;%20x_2%20&#x5C;ge%200%20&#x5C;;%20Integer%20&#x5C;&#x5C;&#x5C;%20x_3%20&#x5C;in%20{0,%201}&#x5C;end{aligned}"/>
  
```julia
  
m4 = Model(Gurobi.Optimizer)
  
@variable(m4, 0 <= x1 <= 10)
@variable(m4, x2 >= 0, Int)
@variable(m4, x3, Bin)
  
@objective(m4, Max, x1 + 2x2 + 5x3)
  
@constraint(m4, constraint1, -x1 + x2 + 3x3 <= -5)
@constraint(m4, constraint2, x1 + 3x2 - 7x3 <= 10)
  
display(m4)
  
optimize!(m4)
  
println("Optimal Solutions:")
for i in index_x
  println("x[$i] = ", value(x[i]))
end
  
println("Dual Variables:")
for j in index_constraints
  println("dual[$j] = ", dual(constraint[j]))
end
  
```
  