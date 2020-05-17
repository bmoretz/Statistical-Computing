# Julia Basics

````julia
using LinearAlgebra;
````




## Vector, Matrix, and Array

````julia
a = [1; 2; 3]
b = [4 5 6]
A = [1 2 3; 4 5 6]

A[1, 3]
A[2, 1]
````


````
4
````



````julia
transpose(A)

A'
````


````
3×2 Adjoint{Int64,Array{Int64,2}}:
 1  4
 2  5
 3  6
````





Column Vectors:

````julia
a = [1; 2; 3]
c = [7; 8; 9]

a'*c

vec
````


````
vec (generic function with 9 methods)
````





Identity

````julia
Matrix(I, 3, 3)
````


````
3×3 Array{Bool,2}:
 1  0  0
 0  1  0
 0  0  1
````



````julia
zeros(4, 1)

zeros(2, 3)

ones(1, 3)
````


````
1×3 Array{Float64,2}:
 1.0  1.0  1.0
````



````julia
B = [1 3 2; 3 2 2; 1 1 1]

inv(B)
````


````
3×3 Array{Float64,2}:
  1.11022e-16   1.0  -2.0
  1.0           1.0  -4.0
 -1.0          -2.0   7.0
````



````julia
B * inv(B)

inv(B)[2, 1]

a = [1; 2; 3]
````


````
3-element Array{Int64,1}:
 1
 2
 3
````



````julia
d = Array{Float64}(undef, 3)
print(d)
````


````
[5.0e-324, 5.0e-324, 0.0]
````



````julia

d[1] = 1
d[2] = 2
d[3] = 3

print(d)
````


````
[1.0, 2.0, 3.0]
````



````julia
p = Array{Float64}(undef, 3, 1)
q = Array{Float64}(undef, 1, 3)

display(p * q)
````


````
3×3 Array{Float64,2}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
````





## Tuples

````julia
pairs = Array{Tuple{Int64, Int64}}(undef, 3)

pairs[1] = (1, 2)
pairs[2] = (2, 3)
pairs[3] = (3, 4)

pairs
````


````
3-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (2, 3)
 (3, 4)
````





Same as:

````julia
pairs = [ (1, 2); (2, 3); (3, 4) ]
````


````
3-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (2, 3)
 (3, 4)
````



````julia
ijk_array = Array{Tuple{Int64, Int64, Int64}}(undef, 3)

ijk_array[1] = (1, 4, 2)
````


````
(1, 4, 2)
````





## Indices and Ranges

````julia
a = [10; 20; 30; 40; 50; 60; 70; 80; 90]

a[1:3]

a[1:3:9]

a[end-2:end]

b = [200; 300; 400]
a[2:4] = b
a

c = collect(1:2:9)
````


````
5-element Array{Int64,1}:
 1
 3
 5
 7
 9
````



````julia
A = [1 2 3; 4 5 6; 7 8 9]

A[:, 2]

A[:, 2:3]

A[3, :]

A[3, :]'

A[3:3, :]

A[2:3, :]
````


````
2×3 Array{Int64,2}:
 4  5  6
 7  8  9
````





## Printing Messages

````julia
println("Hello, World!")
````


````
Hello, World!
````



````julia
a = 123.0

println("The value of a = ", a)
````


````
The value of a = 123.0
````



````julia
println("a is $a, and a-10 is $(a-10).")
````


````
a is 123.0, and a-10 is 113.0.
````



````julia
b = [1; 3; 10]
println("b is $b")
````


````
b is [1, 3, 10]
````



````julia

using Printf

@printf("The %s of a = %f", "value", a)
````


````
The value of a = 123.000000
````



````julia

c = [123.12345    ;
      10.983      ;
       1.9832132  ]

for i in 1:length(c)

      println("c[$i] = $(c[i])")

end
````


````
c[1] = 123.12345
c[2] = 10.983
c[3] = 1.9832132
````



````julia

for i in 1:length(c)

      @printf("c[%d] = %7.3f\n", i, c[i])
end
````


````
c[1] = 123.123
c[2] =  10.983
c[3] =   1.983
````



````julia

str = @sprintf("The %s of a = %f", "value", a)
println(str)
````


````
The value of a = 123.000000
````





## Collections, Dictionary and For-Loop

````julia
for i in 1:5
      println("This is number $i.")
end
````


````
This is number 1.
This is number 2.
This is number 3.
This is number 4.
This is number 5.
````



````julia


for i in 1:5
      if i >= 3
            break
      end

      println("This is number $i.")
end
````


````
This is number 1.
This is number 2.
````



````julia
s = 0
for i in 1:10
      global s += i
end
println(s)
````


````
55
````



````julia
my_keys = ["Zinedine Zidane", "Magic Johnson", "Yuna Kim"]
my_values = ["football", "basketball", "figure skating"]

d = Dict()

for i in 1:length(my_keys)
      d[my_keys[i]] = my_values[i]
end

display(d)
````


````
Dict{Any,Any} with 3 entries:
  "Magic Johnson"   => "basketball"
  "Zinedine Zidane" => "football"
  "Yuna Kim"        => "figure skating"
````



````julia

for (key, value) in d
      println("$key is a $value player.")
end
````


````
Magic Johnson is a basketball player.
Zinedine Zidane is a football player.
Yuna Kim is a figure skating player.
````



````julia

d["Diego Mardona"] = "football"

d
````


````
Dict{Any,Any} with 4 entries:
  "Magic Johnson"   => "basketball"
  "Zinedine Zidane" => "football"
  "Diego Mardona"   => "football"
  "Yuna Kim"        => "figure skating"
````



````julia
links = [ (1, 2), (3, 4), (4, 2) ]
link_costs = [ 5, 13, 8 ]

link_dict = Dict()

for i in 1:length(links)
      link_dict[ links[i] ] = link_costs[ i ]
end

link_dict

for (link, cost) in link_dict
      println("Link $link has cost of $cost.")
end
````


````
Link (1, 2) has cost of 5.
Link (4, 2) has cost of 8.
Link (3, 4) has cost of 13.
````





### Functions

$f(x, y) = 3x + y$

````julia
function f(x, y)
      return 3x + y
end

f(1, 3)
3 * ( f(3, 2) + f(5, 6) )
````


````
96
````





or

````julia
f(x, y) = 3x + y

f(1, 3)
````


````
6
````



````julia
function my_func(n, m)
      a = zeros(n , 1)
      b = ones(m, 1)
      return a, b
end

x, y = my_func(3, 2)

x

y
````


````
2×1 Array{Float64,2}:
 1.0
 1.0
````



````julia
function f(x)
      return x+2
end

function g(x)
      return 3x+3
end
````


````
g (generic function with 1 method)
````



````julia
function f1(x)
      return x + a
end

a = 0

for i in 1:10
      global a = i
      println(f1(1))
end
````


````
2
3
4
5
6
7
8
9
10
11
````



````julia
function f2(x)
      a = 0
      return x+a
end

a = 5
println(f2(1))
````


````
1
````



````julia
println(a)
````


````
5
````



````julia
function f3(x)
      _a = 0
      return x + _a
end

a = 5
println(f3(1))
````


````
1
````



````julia
println(a)
````


````
5
````



````julia
function f4(x, a)
      return x + a
end

a = 5
println(f4(1, a))
````


````
6
````



````julia
println(a)
````


````
5
````





## Random Number Generation

````julia
rand()
````


````
0.7172178308233799
````



````julia
rand(5)
````


````
5-element Array{Float64,1}:
 0.8446117923143865
 0.9035783540039999
 0.04250589086895307
 0.015412071876057087
 0.39324146082594846
````



````julia
rand(4, 3)
````


````
4×3 Array{Float64,2}:
 0.809962  0.507379  0.650565
 0.682745  0.775755  0.0181288
 0.65382   0.906314  0.710279
 0.198878  0.616287  0.491441
````



````julia
rand() * 100
````


````
19.837565637683795
````



````julia
rand(1:10)

randn(2, 3)
````


````
2×3 Array{Float64,2}:
 -0.615127  -3.0847    0.502431
 -0.587144  -0.827289  0.674643
````



````julia
using StatsFuns;

mu = 50; sigma = 3

normpdf(mu, sigma, 52)

normcdf(mu, sigma, 50)

norminvcdf(mu, sigma, 0.5)
````


````
50.0
````





## File I/O

````julia
datafilename = "data.txt"
datafile = open(datafilename)
data = readlines(datafile)
close(datafile)

println(data)
````


````
["This is the first line.", "This is the second line.", "This is the third 
line."]
````



````julia
outputfilename = "results1.txt"
outputfile = open(outputfilename, "w")
print(outputfile, "Majic Johnson")
println(outputfile, " is a basketball player.")
println(outputfile, "Michael Jordan is also a basketball player.")
close(outputfile)
````



````julia
outputfilename = "results2.txt"
outputfile = open(outputfilename, "a")
println(outputfile, "Yuna Kim is a figure skating player.")
close(outputfile)
````



````julia
using CSV; using DataFrames;

csvfilename = "data.csv"
csvdata = CSV.file(csvfilename) |> DataFrame!

csvdata
````


````
5×3 DataFrame
│ Row │ start node │ end node │ link length │
│     │ Int64      │ Int64    │ Float64     │
├─────┼────────────┼──────────┼─────────────┤
│ 1   │ 1          │ 2        │ 2.0         │
│ 2   │ 1          │ 3        │ 4.5         │
│ 3   │ 2          │ 3        │ 6.0         │
│ 4   │ 2          │ 4        │ 3.0         │
│ 5   │ 3          │ 4        │ 5.0         │
````





## Plotting

````julia
using PyPlot
pygui(false)
# Preparing a figure object
fig = figure()

# Data
x = range(0, stop = 2*pi, length = 1000)
y = sin.(3*x)

# Plotting with linewidth and linestyle specified
plot(x, y, color="blue", linewidth=2.0, linestyle="--")

# Labeling the axes
xlabel(L"value of $x$")
ylabel(L"\sin(3x)")

# Title
title("Test plotting")
plt.show()

p = plot(x,y)
xlabel("X")
ylabel("Y")
PyPlot.title("Your Title Goes Here")
grid("on")

p = plot_date(x,y,linestyle="-",marker="None",label="Base Plot") # Basic line plot
````


````
1-element Array{PyObject,1}:
 PyObject <matplotlib.lines.Line2D object at 0x000000005DAE16C8>
````



````julia
using PyPlot

plot([1,2,3,4])

# Data
lower_bound = [4.0, 4.2, 4.4, 4.8, 4.9, 4.95, 4.99, 5.00]
upper_bound = [5.4, 5.3, 5.3, 5.2, 5.2, 5.15, 5.10, 5.05]
iter = 1:8

# Creating a new figure object
fig = figure()

# Plotting two datasets
plot(iter, lower_bound, color="red", linewidth=2.0, linestyle="-",
 marker="o", label=L"Lower Bound $Z^k_L$")
plot(iter, upper_bound, color="blue", linewidth=2.0, linestyle="-.",
 marker="D", label=L"Upper Bound $Z^k_U$")

# Labeling axes
xlabel(L"iteration clock $k$", fontsize="xx-large")
ylabel("objective function value", fontsize="xx-large")

# Putting the legend and determining the location
legend(loc="upper right", fontsize="x-large")

# Add grid lines
grid(color="#DDDDDD", linestyle="-", linewidth=1.0)
tick_params(axis="both", which="major", labelsize="x-large")

# Title
title("Lower and Upper Bounds")

# Save the figure as PNG and PDF
savefig("plot2.png")
savefig("plot2.pdf")

# Closing the figure object
close(fig)
````


