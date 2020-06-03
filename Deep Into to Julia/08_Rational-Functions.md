---
title: 'Polynomial Roots'
author: 'Brandon Moretz'
date: '1st June 2020'
---  
  
  
#  Rational Functions
  
  
```julia
using CalculusWithJulia   # to load the `Plots` and `SymPy` packages
  
f(x) = (x-1)^2 * (x-2) / ((x+3)*(x-3) )
plot(f, -10, 10)
```
  
```julia
plot(f, -100, 100)
```
  
```julia
@vars x real=true
a = (x-1)^2 * (x-2)
b = (x+3)*(x-3)
```
  
```julia
q, r = divrem(a, b)
```
  
```julia
q + r/b
```
  
```julia
f(x) = (x-1)^2 * (x-2) / ((x+3)*(x-3))  # as a function
p = f(x)                                # a symbolic expression
apart(p)
```
  
```julia
plot(apart(p) - (x - 4), 10, 100)
```
  
```julia
cancel(p)
```
  