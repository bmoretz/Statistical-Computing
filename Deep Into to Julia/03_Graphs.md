# Graphs

````julia
using CalculusWithJulia
````





$f(x) = 1 - \frac{x^2}{2}$

````julia
f(x) = 1 - x^2/2
plot(f, -3, 3)
````


![](figures/03_Graphs_2_1.png)

````julia
plot(sin, 0, 2pi)
````


![](figures/03_Graphs_3_1.png)



$f(x) = (1 + x^2)^{-1}$

````julia
f(x) = 1 / (1 + x^2)
plot(f, -3, 3)
````


![](figures/03_Graphs_4_1.png)



$f(x) = e^{-x^2/2}$

````julia
f(x) = exp(-x^2/2)
plot(f, -2, 2)
````


![](figures/03_Graphs_5_1.png)

````julia
plot(x -> cos(x), 0, pi/2)
````


![](figures/03_Graphs_6_1.png)

````julia
f(x) = tan(x)
plot(f, -10, 10)
````


![](figures/03_Graphs_7_1.png)

````julia
f(x) = sin(x)
xs = range(0, 2pi, length = 10)
ys = f.(xs)

plot(xs, ys)
````


![](figures/03_Graphs_8_1.png)

````julia
f(x) = 1/x
xs = range(-1, 1, length = 251)
ys = f.(xs)
ys[xs .== 0.0] .= NaN

plot(xs, ys)
````


![](figures/03_Graphs_9_1.png)

````julia
g(x) = abs(x) < .05 ? NaN : f(x)
plot(g, -1, 1)
````


![](figures/03_Graphs_10_1.png)

````julia
xs = range(-pi, pi, length = 100)
ys = sin.(xs)

plot(xs, ys)
````


![](figures/03_Graphs_11_1.png)

````julia
plot(-xs, ys)
````


![](figures/03_Graphs_12_1.png)

````julia
plot(ys, xs)
````


![](figures/03_Graphs_13_1.png)

````julia
xs = range(-pi/2, pi/2, length = 100)
ys = [sin(x) for x in xs]

plot(ys, xs)
````


![](figures/03_Graphs_14_1.png)

````julia
f(x) = 1 - x^2/2
plot([cos, f], -pi/2, pi/2)
````


![](figures/03_Graphs_15_1.png)

````julia
f(x) = x^5 - x + 1
plot([f, zero], -1.5, 1.4)
````


![](figures/03_Graphs_16_1.png)

````julia
plot(f, -1.5, 1.4)
plot!(zero)
````


![](figures/03_Graphs_17_1.png)

````julia
f(x) = x*(x-1)
plot(f, -1, 2)
scatter!([0, 1], [0,0])
````


![](figures/03_Graphs_18_1.png)

````julia
plot(f, title="plot of x*(x-1)",
    xlab = "x axis", ylab = "y axis")

plot(f, linewidth = 5)

plot(f, legend=false)

plot(f, linestyle=:dash)
plot(f, linestyle=:dot)
plot(f, linestyle=:dashdot)

scatter(f, marker = :square, legend=false)
````


![](figures/03_Graphs_19_1.png)



## Parametric Graphs

````julia
f(x) = cos(x); g(x) = sin(x)
xs = range(0, 2pi, length = 100)
plot(f.(xs), g.(xs))
````


![](figures/03_Graphs_20_1.png)

````julia
plot(f, g, 0, 2pi)
````


![](figures/03_Graphs_21_1.png)

````julia
g(x) = x^2
f(x) = x^3
plot([g, f], 0, 25)
````


![](figures/03_Graphs_22_1.png)

````julia
xs = range(0, 5, length = 100)
plot(g, f, 0, 25)
plot(f, g, 0, 25)
````


![](figures/03_Graphs_23_1.png)

````julia
g(x) = x - sin(x)
f(x) = x^3
plot(g, f, -pi/2, pi/2)
````


![](figures/03_Graphs_24_1.png)

````julia
g(x) = x - sin(x)
f(x) = x^3

plot(g)
plot!(f)

plot(g, f, -pi/2, pi/2)
````


![](figures/03_Graphs_25_1.png)

````julia
R, r, rho = 1, 1/4, 1/4

g(t) = (R-r) * cos(t) + rho * cos((R-r)/r * t)
f(t) = (R-r) * sin(t) - rho * sin((R-r)/r * t)

plot(g, f, 0, max((R-r)/r, r/(R-r))*2pi)
````


![](figures/03_Graphs_26_1.png)

````julia
f(x) = x^3 - x
plot([f, zero], -2, 2)
````


![](figures/03_Graphs_27_1.png)

````julia
f(x) = x^3 - x

plot([f, zero])
````


![](figures/03_Graphs_28_1.png)



Given,

$f(x) = 3x^4 + 8x^3 - 18x^2$

Find the point at which f(x) is the smallest.

````julia
f(x) = 3x^4 + 8x^3 - 18x^2

xs = range(-4, -2, length = 100)
ys = f.(xs)

plot(f)

xs[argmin(ys)]
````


````
-2.98989898989899
````





$f(x) 3x^4 + 8x^3 - 18x^2$

When is it increasing?

````julia
f(x) = 3x^4 + 8x^3 - 18x^2

plot(f)
````


![](figures/03_Graphs_30_1.png)



$f(x) = \frac{(x^3 - 2x)}{2x^2 - 10}$

is a rational function with issues When

$2x^2 = 10$ or $x = ±\sqrt{5}$

````julia
f(x) = (x^3 - 2x)/(2x^2 - 10)
plot([f, zero], -5, 5)
````


![](figures/03_Graphs_31_1.png)

````julia
f(x) = x <= 10 ? 35.0 : 35.0 + 4.0 * (x-10)

plot(f, 0, 20)
hline!([55])
````


![](figures/03_Graphs_32_1.png)

````julia
f(x) = cos(x); g(x) = x

plot([f, g])
vline!([.75])
````


![](figures/03_Graphs_33_1.png)

````julia
f(x) = log(x)-2
plot([f, zero],0, 10)
vline!([7.5])
````


![](figures/03_Graphs_34_1.png)

````julia
xs = range(0, 1, length=250)
f(x) = sin(500*pi*x)
plot(xs, f.(xs))

plot(f, 0, 1)
````


![](figures/03_Graphs_35_1.png)

````julia
function trimplot(f, a, b, c=20; kwargs...)
   fn = x -> abs(f(x)) < c ? f(x) : NaN
   plot(fn, a, b; kwargs...)
end

f(x) = 1/x
plot(f, -1, 1)
trimplot(f, -1, 1)
````


![](figures/03_Graphs_36_1.png)

````julia
R, r, rho = 1, 3/4, 1/4

f(t) = (R-r) * cos(t) + rho * cos((R-r)/r *t)
g(t) = (R-r) * sin(t) + rho * sin((R-r)/r *t)

plot(f, g, 0, max((R-r)/r, r/(R-r))*2pi, aspect_ratio=:equal)
````


![](figures/03_Graphs_37_1.png)

````julia
function spirograph(R, r, rho)
  f(t) = (R-r) * cos(t) + rho * cos((R-r)/r * t)
  g(t) = (R-r) * sin(t) - rho * sin((R-r)/r * t)

  plot(f, g, 0, max((R-r)/r, r/(R-r))*2pi, aspect_ratio=:equal)
end

spirograph(1, 3/4, 1/4)
````


![](figures/03_Graphs_38_1.png)

````julia
spirograph(1, 1/2, 1/4)
````


![](figures/03_Graphs_39_1.png)

````julia
spirograph(1, 1/4, 1)
````


![](figures/03_Graphs_40_1.png)

````julia
spirograph(1, 1/8, 1/4)
````


![](figures/03_Graphs_41_1.png)
