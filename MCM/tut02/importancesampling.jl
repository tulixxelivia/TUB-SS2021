using Distributions
using LinearAlgebra
using SpecialFunctions
using Random
using Plots

"""
    importancesampling(target, proposal, n)

## Arguments

- `target` : Distribution to target
- `proposal` : proposal distribution
- `n` : number of samples

"""
function importancesampling(target, proposal, n)
    x = rand(proposal, n)
    return x, exp.(logpdf(target, x) - logpdf(proposal, x)) # returns the samples and weights
end

function is_expectation(f, target, proposal, n)
    x, w = importancesampling(target, proposal, n)
    return sum(w .* f.(x)) / n
end

function naive_expectation(f, p, n)
    x = rand(p, n)
    return sum(f.(x)) / n
end
## Setup 
a = 5
b = 4
f(x) = x^(1 - a)
p = Gamma(a, b)
params = (m = (a - 1) * b, v = sqrt(a * b^2)) # Good
params = (m = a * b, v = sqrt(a * b^2)) # Bad
# params = (m = 0, v = sqrt(a * b^2)) # Better
# params = (m = 2 * a * b, v = sqrt(a * b^2)) # Worse
# q = Normal(params...)
q = Exponential(1 / b)

true_val = b^(1.0 - a) / gamma(a) #  true value of the expectation
 
# Plotting
N = 1000

plt1 = histogram(rand(q, N), normalize=true, label="proposal samples", alpha=0.5)
histogram!(plt1, rand(p, N), normalize=true, label="target samples", alpha=0.5)
plot!(plt1, x->pdf(q, x); label="q(x)", lw=3.0, color=1)
plot!(plt1, x->pdf(p, x); label="p(x)", lw=3.0, color=2)
plot!(plt1, x->f(x); label="f(x)", lw=3.0, color=:black, ylims=(0, ylims(plt1)[2]), xlims=xlims(plt1))

# Estimator
N = 1000
M = 100
naive_estimations = [naive_expectation(f, p, N) for i in 1:M]
is_estimations = [is_expectation(f, p, q, N) for i in 1:M]
plt2 = vline([true_val], label="True estimate", lw=5.0, color=:black)
histogram!(plt2, naive_estimations, normalize=true, label="Naive expectation", alpha=0.3, xlims=(0, true_val * 2), color=2)
histogram!(plt2, is_estimations, normalize=true, label="IS expectation", alpha=0.3, color=1)

plot(plt1, plt2) |> display

# Inspect stats about estimators
@show mean(naive_estimations), var(naive_estimations)
@show mean(is_estimations), var(is_estimations)
@show true_val


## A different problem now!

a = 5
b = 10
m = (a + b) / 2
f(x) = exp(-abs(x-m)/10)
f(x) = a < x < b
p = Normal(0, 1)
q = Normal(m, 1)


plt1 = histogram(rand(q, N), normalize=true, label="proposal samples", alpha=0.5)
histogram!(plt1, rand(p, N), normalize=true, label="target samples", alpha=0.5)
plot!(plt1, x->pdf(q, x); label="q(x)", lw=3.0, color=1)
plot!(plt1, x->pdf(p, x); label="p(x)", lw=3.0, color=2)
plot!(plt1, x->f(x); label="f(x)", lw=3.0, color=:black, ylims=(0, 1.1), xlims=xlims(plt1))

true_val = cdf(p, b) - cdf(p, a)
N = 1000
M = 100
is_expectation(f, p, q, N)
naive_expectation(f, p, N)
naive_estimations = [naive_expectation(f, p, N) for i in 1:M]
is_estimations = [is_expectation(f, p, q, N) for i in 1:M]
plt2 = vline([true_val], label="True estimate", lw=5.0, color=:black)
histogram!(plt2, naive_estimations, normalize=true, label="Naive expectation", alpha=0.3, color=2)
histogram!(plt2, is_estimations, normalize=true, label="IS expectation", alpha=0.3, color=1)

@show mean(naive_estimations), var(naive_estimations)
@show mean(is_estimations), var(is_estimations)
@show true_val

plot(plt1, plt2) |> display
