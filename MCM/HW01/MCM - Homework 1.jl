### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ cb89d54f-7a8f-4f97-9405-899971e70483
begin
	using Pkg; Pkg.add(["Distributions", "LinearAlgebra", "Plots", "PlutoUI"])
	using Distributions
	using LinearAlgebra
	using Plots
	using PlutoUI
	default(legendfontsize = 15.0, linewidth = 2.0)
end

# ╔═╡ c9b39b6b-7489-428d-b3fc-d4c5f1f180c4
# TableOfContents()

# ╔═╡ d7a6e8e6-9c69-11eb-1c3b-1d3fdd4450ec
md"""# Problem Sheet 1"""

# ╔═╡ 7e4c9518-2eb8-426b-a5e6-c9d4b2bd6499
md"""## Inverse Transformation Method: Cauchy Distribution"""

# ╔═╡ b8a8b63b-77d1-4e5d-bd00-a5b9f0e34925
md"""The probability density function (pdf) of a cauchy distribution is defined
```math
p(x|x_0, \gamma) = \frac{1}{\pi}\left(\frac{\gamma}{(x-x_0)^2+\gamma^2}\right)
```
 
where $x_0$ is the location of the mode and $\gamma$ is a shape parameter.

Use the inverse transformation method to find a function $f(u|x_0, \gamma)$ which generates cauchy distributed samples if $u$ is uniformly distributed over $[0,1]$."""

# ╔═╡ 2f8775b1-9d69-4b9a-82a8-2429768604b2
md"""
*Fill in your answer here or on paper!*
"""

# ╔═╡ d38c81ac-4af0-4ead-991a-f94bfccc11df
md"""
γ = $(@bind γ Slider(0.01:0.01:5; default=2.0, show_value=true))

x₀ = $(@bind x₀ Slider(-5:0.01:5; default=3.0, show_value=true))
"""

# ╔═╡ 33ade8d2-898d-481d-a1ac-b7fd393b239f
true_d = Cauchy(x₀, γ) # https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Cauchy

# ╔═╡ 6ba89d7c-2ad8-4f26-b569-2aa4fe86cb56
md"We sample N uniform variables and compute the inverse transform on the uniform variables"

# ╔═╡ 40285555-ffad-4fe8-9f9f-2ea9802b9890
function f(u, x₀, γ)
	## !! CODE MISSING !! ##
	## Fill in the function correctly transforming your uniform random variables here
	## !! CODE MISSING !! ##
end;

# ╔═╡ 53d4924c-2435-49f3-a1e8-1b41383492a2
begin
	N = 1000
	u = rand(N)
	x = f.(u, x₀, γ) # We broadcast the transformation on all u
end;

# ╔═╡ 7d80aa2e-b447-44a8-8d58-18a1f6692c43
true_x = rand(true_d, N); # We also sample from the method in Distributions.jl

# ╔═╡ 5a643a9d-6064-446e-bda0-2701bbefb7eb
## Plotting some histograms
begin
	bins = range(-10 + x₀, stop = 10 + x₀,length = 30) # We create a range for the bins (adapted to x₀)

# We plot the histogram of the samples, normalize it and compare it to the true pdf
	if x[1] !== nothing
		h = histogram(x, label="", title = "Inverse Transform", bins=bins, normalize = true, ylims = (0, 0.2))
		plot!(h, x->pdf(true_d, x), linewidth = 3.0, label = "")
		true_h = histogram(true_x, label = "", title = "True Distribution", bins = bins, normalize = true, ylims = (0, 0.2))
		plot!(true_h, x->pdf(true_d, x), linewidth = 3.0, label = "")
		plot(h, true_h)
	end
end

# ╔═╡ 0a881e39-fb8e-4853-bcf0-58bd172d34f0
md"## Polar Box-Muller"

# ╔═╡ eb89f195-759c-4423-bb1d-61e44eb4a2f3
md"""
A computational more efficient version of the Box-Muller transformation makes use of random numbers $z_1$, $z_2$ which are uniformly distributed in the unit circle.

We can generate these by drawing $z_1$ and $z_2$ from the uniform distribution over $[-1,1]$ and rejecting the pairs until $z_1^2+z_2^2 \leq 1$ is true.


- **Show** that
```math
\begin{align}
 y_1 &= z_1\sqrt{\frac{-2\ln(r^2)}{r^2}}~(1)\\
 y_2 &= z_2\sqrt{\frac{-2\ln(r^2)}{r^2}}~(2),
\end{align}
```
with $r^2=z_1^2 + z_2^2$ ,have the joint distribution
```math
\begin{equation}
 p(y_1,y_2) = \left(\frac{1}{\sqrt{2\pi}}\exp\left(\frac{-y_1^2}{2}\right)\right)\left(\frac{1}{\sqrt{2\pi}}\exp\left(\frac{-y_2^2}{2}\right)\right)
\end{equation}
```
and therefore each has a standard normal distributio (Gaussian distribution with zero mean and unit variance).

!!! tip
	Use the two first equations (1) and (2) to get a formula for $r^2$ depending only on $y_1$ and $y_2$. Then you can easily get the inverted function $z_{1/2}(y_1,y_2)$.

	You might want to use *Mathematica* or *Maple* for computing some of the derivates. If not you can use
	```math
	\begin{equation}
	 [y_1^4-2y_2^2+y_1^2y_2^2][y_2^4-2y_1^2+y_1^2y_2^2] - [y_1y_2(2+y_1^2+y_2^2)]^2 = - 2(y_1^2+y_2^2)^3
	\end{equation}
	```
"""

# ╔═╡ 5495ebee-8563-452b-85d4-1cbd21e54bf6
md"""
Fill in your answers here or on paper
"""

# ╔═╡ f3eb69ca-e06e-4399-8ffa-60cc998ce3d7
md"## Rejection Sampler: General"

# ╔═╡ 6679410c-9952-4815-9098-436ed3393ea5
md"""
We want to use a rejection sampler to sample from a target distribution $p(x)$. As a proposal distribution we can choose between:

- (1) $q_1(x)$ , with $c_1 = 1.5$, 
- (2) $q_2(x)$, with $c_2 = 2.0$,
- (3) $q_3(x)$, with $c_3 = 4.0$,
 
where for $i=1,2,3$ the constant $c_i$ is the smallest number which fulfills $c_i q_i(x) \geq p(x) ~\forall x\in \mathcal{R}$.

We know that on average it takes $6 \cdot 10^{-4}$, $4 \cdot 10^{-4}$ and $3 \cdot 10^{-4}$ seconds to get one sample from $q_1$, $q_2$ and $q_3$, respectively.

  With this information, **which proposal distribution should you use and why?**
"""

# ╔═╡ f0257709-48d1-4514-bb0e-52694ea404c7
md"""
*Fill in you answer here or on paper*
"""

# ╔═╡ 0b18b432-4f5b-4696-9d83-ca776b8a2a5f
md"## Rejection Sampler: Rayleigh distribution"

# ╔═╡ 60354c6b-eddc-46ad-8d50-b54db6f4ca76
md"""
We want to use a rejection sampler to sample from a Rayleigh distribution. The pdf of a Rayleigh distribution is
```math
\begin{equation}
     p(x|\sigma_R) = 
     \begin{cases}
      \frac{x}{\sigma_R^2}\exp\left(-\frac{x^2}{2\sigma_R^2}\right), &x \geq 0\\
      0, & x < 0. 
     \end{cases}
\end{equation}
```
We want to use a Gaussian distribution as a proposal distribution and set its mean to $\sigma_R$ (which is the mode of the Rayleigh distribution) and its standard deviation to $\sigma_G$.
"""

# ╔═╡ f2065ee1-8137-476f-bc0a-62e57ac9add2
md"""
### a) Show that if $\sigma_G \geq \sigma_R$, $c$ has to be at least

```math
\begin{equation}
     \sqrt{\frac{2\pi}{\exp(1)}}\frac{\sigma_G}{\sigma_R}
\end{equation}
```
"""

# ╔═╡ df422247-32fd-43e8-94cb-45e684103271
md"""
*Fill in your answer here or on paper*
"""

# ╔═╡ dc2dc437-592a-438d-831c-bfc8c7b3b527
md"""
For $\sigma_G \geq \sigma_R$ only $x_1$ is positive.
```math
\begin{align*}
		\frac{d^2c}{dx^2}(x) &=\left[
		\left(\frac{2x(1-\frac{\sigma_G^2}{\sigma_R^2})-\sigma_R}{\sigma_G\sigma_R^2}\right)
		+\left(\frac{2x(\sigma_R^2-\sigma_G^2)-2\sigma_R^3}{2\sigma_G^2\sigma_R^2}\right)\left(\frac{x^2(1-\frac{\sigma_G^2}{\sigma_R^2})-x\sigma_R + \sigma_G^2}{\sigma_G\sigma_R^2}\right)
		\right]\\
		&\times \sqrt{2\pi}\exp\left(\frac{x^2(\sigma_R^2-\sigma_G^2) - 2x\sigma_R^3+\sigma_R^4}{2\sigma_G^2\sigma_R^2}\right)\\
		\end{align*}
```
```math
\begin{align*}
		\text{sign}\left(\frac{d^2c}{dx^2}(\sigma_R)\right) &=\text{sign}(
		\left(\frac{2\sigma_R(1-\frac{\sigma_G^2}{\sigma_R^2})-\sigma_R}{\sigma_G\sigma_R^2}\right)
		+\left(\frac{2\sigma_R(\sigma_R^2-\sigma_G^2)-2\sigma_R^3}{2\sigma_G^2\sigma_R^2}\right)\left(\frac{\sigma_R^2(1-\frac{\sigma_G^2}{\sigma_R^2})-\sigma_R^2 + \sigma_G^2}{\sigma_G\sigma_R^2}\right))\\
		&= \text{sign}\left(
		\left(\frac{\sigma_R-\frac{\sigma_G^2}{\sigma_R}}{\sigma_G\sigma_R^2}\right)
		+\left(\frac{-2\sigma_G^2\sigma_R)}{2\sigma_G^2\sigma_R^2}\right)\left(\frac{\sigma_R^2-\sigma_G^2-\sigma_R^2 + \sigma_G^2}{\sigma_G\sigma_R^2}\right)\right)\\
		&= \text{sign}
		\left(\frac{\sigma_R^2-\sigma_G^2}{\sigma_G\sigma_R^3}\right)\\
		&= -1, \text{for }\sigma_G \geq \sigma_R.
\end{align*}
```
"""

# ╔═╡ a5616e95-9476-45db-9773-edbe3ec2639d
md"""
σr = $(@bind σr Slider(0.01:0.1:5; default=1.0, show_value=true))
σg = $(@bind σg Slider(0.01:0.1:5; default=2.0, show_value=true))
"""

# ╔═╡ ec65ffcd-ebc4-414a-9ae3-90f6640472f0
min_c = 5.0 ## !! CODE MISSING !! ## Write here what hte minimum value of c should be

# ╔═╡ 387b51cb-456a-4ae5-bfed-f1390ade3503
md"c = $(@bind c Slider(0.01:0.1:10; default=min_c, show_value=true))"

# ╔═╡ 7875a5b7-c24b-455a-99e4-c77799b41638
begin
	xrange = range(-1, 8, length = 200)
	d_rayleigh = Rayleigh(σr) # We create a Rayleigh distribution
	plot(xrange, pdf.(d_rayleigh, xrange); label = "Rayleigh", xlabel="x", ylabel="p(x)") # We plot the pdf
	d_gauss = Normal(σr, σg)
	plot!(xrange, pdf.(d_gauss, xrange), label = "Gaussian")
	plot!(xrange, c * pdf.(d_gauss, xrange), label = "c * Gaussian")
end

# ╔═╡ f010c1a2-70cd-4aa2-bd07-d1e444f13196
md"Let's plot the same thing but in log-space"

# ╔═╡ c95311f6-0a0d-4f62-aaf9-e772f813b87f
begin
	logxrange = range(0, 100, length = 100)
	plot(logxrange, logpdf.(d_rayleigh, logxrange), label = "Rayleigh", xlabel="x", ylabel="log p(x)")
	plot!(logxrange, logpdf.(d_gauss, logxrange), label = "Gaussian")
	plot!(logxrange, log(c) .+ logpdf.(d_gauss, logxrange), label = "c * Gaussian", style = :dot)
end

# ╔═╡ 5db06db1-b8f7-498d-84ac-acc068d9f1e8
md"""
### (b) It is obvious that $c$ is minimal for $\sigma_G = \sigma_R$. Why is it not possible to use a Gaussian distribution with $\sigma_G < \sigma_R$ as a proposal?
"""

# ╔═╡ f00c54a3-6d25-4ce0-9bcd-771075a35ed3
md"""
*Fill in your answer here or on paper*
"""

# ╔═╡ 93714a72-52dc-4fab-b8ce-76e5cd660b5a
md"""
### (c) Program a rejection sampler for a Rayleigh distribution with an arbitary $\sigma_R$ and a Gaussian distribution with standard deviation and mean $\sigma_R$ and the minimal $c$ as computed in (a). Plot the acceptance rate against the number of drawn samples to show that it converges to $\frac{1}{c}$.
"""

# ╔═╡ d859df4a-433f-4422-8737-fd42b8cab5fb
md"""
#### Solution
"""

# ╔═╡ b2262967-ce94-439c-a049-0e8eec66dbff
begin
	N_samples = 5000
	accepted = falses(N_samples)
	acceptance_rate = zeros(N_samples)
	samples = zeros(N_samples)
	for i in 1:N_samples
		x = rand(d_gauss)
		## !! CODE MISSING !! 
		## accepted[i] = ? # Fill in the rejection sampler here, accepted[i] takes a boolean indicating if the sample is accepted
		## !! CODE MISSING !! 
		if accepted[i]
			samples[i] = x
		end
		acceptance_rate[i] = sum(accepted[1:i]) / i
	end
end

# ╔═╡ d274fb70-05dd-413c-a50e-f1a3b3b72eb8
begin
	plot(acceptance_rate, label = "Acceptance rate", lw = 3.0, ylims = (0,1))
	plot!(1 / c * ones(N_samples), label = "1/c", lw = 3.0)
end

# ╔═╡ e5b6a256-f70c-401d-b2e6-81bf2f369622
begin
	if acceptance_rate[end] > 0
		histogram(samples[samples .!= 0]; normalized=true, label="")
		plot!(x->pdf(d_rayleigh, x); lw=4.0, label="")
	end
end

# ╔═╡ Cell order:
# ╠═cb89d54f-7a8f-4f97-9405-899971e70483
# ╠═c9b39b6b-7489-428d-b3fc-d4c5f1f180c4
# ╟─d7a6e8e6-9c69-11eb-1c3b-1d3fdd4450ec
# ╟─7e4c9518-2eb8-426b-a5e6-c9d4b2bd6499
# ╟─b8a8b63b-77d1-4e5d-bd00-a5b9f0e34925
# ╟─2f8775b1-9d69-4b9a-82a8-2429768604b2
# ╟─d38c81ac-4af0-4ead-991a-f94bfccc11df
# ╟─33ade8d2-898d-481d-a1ac-b7fd393b239f
# ╟─6ba89d7c-2ad8-4f26-b569-2aa4fe86cb56
# ╠═40285555-ffad-4fe8-9f9f-2ea9802b9890
# ╠═53d4924c-2435-49f3-a1e8-1b41383492a2
# ╠═7d80aa2e-b447-44a8-8d58-18a1f6692c43
# ╟─5a643a9d-6064-446e-bda0-2701bbefb7eb
# ╟─0a881e39-fb8e-4853-bcf0-58bd172d34f0
# ╟─eb89f195-759c-4423-bb1d-61e44eb4a2f3
# ╟─5495ebee-8563-452b-85d4-1cbd21e54bf6
# ╟─f3eb69ca-e06e-4399-8ffa-60cc998ce3d7
# ╟─6679410c-9952-4815-9098-436ed3393ea5
# ╟─f0257709-48d1-4514-bb0e-52694ea404c7
# ╟─0b18b432-4f5b-4696-9d83-ca776b8a2a5f
# ╟─60354c6b-eddc-46ad-8d50-b54db6f4ca76
# ╟─f2065ee1-8137-476f-bc0a-62e57ac9add2
# ╟─df422247-32fd-43e8-94cb-45e684103271
# ╟─dc2dc437-592a-438d-831c-bfc8c7b3b527
# ╟─a5616e95-9476-45db-9773-edbe3ec2639d
# ╠═387b51cb-456a-4ae5-bfed-f1390ade3503
# ╠═ec65ffcd-ebc4-414a-9ae3-90f6640472f0
# ╟─7875a5b7-c24b-455a-99e4-c77799b41638
# ╟─f010c1a2-70cd-4aa2-bd07-d1e444f13196
# ╟─c95311f6-0a0d-4f62-aaf9-e772f813b87f
# ╟─5db06db1-b8f7-498d-84ac-acc068d9f1e8
# ╟─f00c54a3-6d25-4ce0-9bcd-771075a35ed3
# ╟─93714a72-52dc-4fab-b8ce-76e5cd660b5a
# ╟─d859df4a-433f-4422-8737-fd42b8cab5fb
# ╟─b2262967-ce94-439c-a049-0e8eec66dbff
# ╟─d274fb70-05dd-413c-a50e-f1a3b3b72eb8
# ╟─e5b6a256-f70c-401d-b2e6-81bf2f369622
