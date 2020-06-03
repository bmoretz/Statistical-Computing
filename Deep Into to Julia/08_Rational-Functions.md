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
  
```julia
p = (x^5 - 2x^4 + 3x^3 - 4x^2 + 5) / (5x^4 + 4x^3 + 3x^2 + 2x + 1)
apart(p)
```
  
```julia
a = 5x^3 + 6x^2 + 2
b = x - 1
q, r = divrem(a, b)
```
  
```julia
plot(a/b, -3, 3)
plot!(q)
```
  
```julia
plot(a/b, 5, 10)
plot!(q)
```
  
```julia
x = symbols("x")
p = (x-1)*(x-2)
q = (x-3)^3 * (x^2-x-1)
apart(p/q)
```
  
```julia
plot(1/x, -1, 1)
```
  
```julia
f(x) = (x-1)^2 * (x-2) / ((x+3)*(x-3))
f(3), f(-3)
```
  
```julia
f(x) = (x-1)^2 * (x-2) / ((x+3)*(x-3) )
plot(f, -2.9, 2.9)
```
  
```julia
plot(f, -5, 5, ylims=(-20, 20))
```
  
```julia
function trim_plot(f, a, b, c=20; kwargs...)
   fn = x -> abs(f(x)) < c ? f(x) : NaN
   plot(fn, a, b; kwargs...)
end
```
  
```julia
trimplot(f, -25, 25, 30)
```
  
```julia
function sign_chart(f, a, b)
   xs = range(a, stop=b, length=200)
   ys = f.(xs)
   cols = [fx < 0 ? :red  : :blue for fx in ys]
   plot(xs, ys, color=cols, linewidth=5, legend=false)
   plot!(zero)
   end
```
  
```julia
f(x) = x^3 - x
signchart(f, -3/2, 3/2)
```
  
```julia
sin_p(x) = (x - (7/60)*x^3) / (1 + (1/20)*x^2)
tan_p(x) = (x - (1/15)*x^3) / (1 - (2/5)*x^2)
plot(sin, -pi, pi)
plot!(sin_p)
```
  
```julia
plot(tan, -pi/2 + 0.2, pi/2 - 0.2)
plot!(tan_p)
```
  
```julia
t = symbols("t")
fr(t) = 50t^2 / (t^3 + 20)
  
fr(1)
plot(fr, 0, 60)
```
  
```julia
fr(24)
```
  
```julia
plot(fr, 0, 60)
vline!([3.5])
```
  