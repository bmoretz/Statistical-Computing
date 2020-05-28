# Polynomial Roots

````julia
using CalculusWithJulia;
````



````julia
f(x) = x^2 - x
plot(f, -2, 2)
plot!(zero)
````


![](figures/07_PolynomialRoots_2_1.png)

````julia
f(x) = log(x)
plot(f, 0, 100)
plot!(zero)

plot([0, 20], [25, 25], color =:green)
````


![](figures/07_PolynomialRoots_3_1.png)
