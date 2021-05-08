using Distributions, Random
using Plots; pyplot()
using ProgressMeter

inverselaplacecdf(u, m, beta) = m - beta * sign.(u) * log.(1.0 .- 2.0 * abs.(u))

"""
    sampleLaplace(alpha,xm,n,drawInterval)

Example of the inverse transformation method (Laplace function distribution)
- mu mean parameter
- beta  shape parameter
- n number of samples
- drawInterval number of samples after which the histogramms are drawn
"""
function sampleLaplace(mu = 0, beta = 1, n = 1000; drawInterval = 50)

    binsuniform = range(-0.5, 0.5, length = 20)
    binslaplace = range(mu - 3 * 2 * beta ^2, mu + 3 * 2 * beta^2, length = 100)
    uniformSamples = zeros(n)
    laplaceSamples = zeros(n)
    d_uniform = Uniform(-0.5, 0.5)
    d_laplace = Laplace(mu, beta)
    anim = Animation()
    @showprogress for i = 1:n
        newUniformSample = rand(d_uniform)
        # The inverse of the cdf:
        newLaplaceSample = inverselaplacecdf(newUniformSample, mu, beta)
        uniformSamples[i] = newUniformSample
        laplaceSamples[i] = newLaplaceSample
        if i % min(drawInterval,10^(floor(log10(i)))) == 0
            p1 = histogram(uniformSamples[1:i]; bins=binsuniform, normalize=true, label="", lw=0.0, title="N = $i")
            plot!(x->pdf(d_uniform, x); label="U(-0.5,0,5)", lw = 3.0)
            p2 = histogram(laplaceSamples[1:i]; bins=binslaplace, normalize=true, label="", lw=0.0)
            plot!(x->pdf(d_laplace, x); label="Laplace", lw=3.0)
            plot(p1, p2)
            frame(anim)
        end
    end
    return gif(anim; fps=5)
end


sampleLaplace()
