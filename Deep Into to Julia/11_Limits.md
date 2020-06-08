---
title: Limits
author: 'Brandon Moretz'
date: '8st June 2020'
---  
  
  
#  Limits
  
  
```julia
using CalculusWithJulia  # to load Plots and SymPy
```
  
```julia
x = 1/10000
(1 + x)^(1/x)
```
  
```julia
f(x) = (1 + x)^(1/x)
xs = [1/10^i for i in 1:10]
[xs f.(xs)]
```
  
```julia
x = 1
sin(x) / x
```
  
```julia
f(x) = sin(x)/x
plot(f, -pi/2, pi/2, legend=false)
```
  
```julia
plot([sin, x -> x], -pi/2, pi/2)
```
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;lim_{x→2}&#x5C;frac{x^2%20-%205x%20+%206}{x^2%20+%20x%20-%206}"/>
  
```julia
f(x) = (x^2 - 5x + 6) / (x^2 + x - 6)
c = 2
f(c)
```
  
```julia
c, delta = 2, 1
plot(f, c - delta, c + delta)
```
  
```julia
f(x) = x == 2.0 ? -0.2 : (x^2 - 5x + 6) / (x^2 + x - 6)
```
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;lim_{x→25}&#x5C;frac{&#x5C;sqrt{x}%20-%205}{&#x5C;sqrt{x%20-%2016}%20-%203}"/>
  
```julia
f(x) = (sqrt(x) - 5)/(sqrt(x - 16) - 3)
c = 25
f(c)
  
plot(f, 0, 40)
```
  
```julia
hs = [1/10^i for i in 1:8]
  
xs = c .+ hs
ys = f.(xs)
```
  
```julia
[xs ys]
```
  
```julia
xs = c .- hs
ys = f.(xs)
[xs ys]
```
  
```julia
c = 1
f(x) = x^x
ys = [(f(c + h) - f(c))/h for h in hs]
[hs ys]
```
  
```julia
ys = [(f(c + h) - f(c))/h for h in -hs]
[-hs ys]
```
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;lim_{x→0}&#x5C;frac{1%20-%20cos(x)}{x^2}"/>
  
```julia
f(x) = (1 - cos(x))/x^2
f(0)
  
c = 0
xs = c .+ hs
ys = [f(x) for x in xs]
[xs ys]
```
  
```julia
y1s = [1 - cos(x) for x in xs]
y2s = [x^2 for x in xs]
[xs y1s y2s]
```
  
```julia
@vars x real=true
f(x) = (1 - cos(x))/x^2
limit(f(x), x => 0)
```
  
```julia
limit(f, 0)
```
  
```julia
limit( (2sin(x) - sin(2x)) / (x - sin(x)), x => 0)
```
  
```julia
f(x) = (exp(x) - 1 - x)/x^2
  
limit(f, 0)
```
  