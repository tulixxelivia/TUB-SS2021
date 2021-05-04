using Distributions, Random
using Plots; pyplot()

inverseparetocdf(u, alpha, xm) = xm / (u ^ (1 / alpha))

"""
    samplePareto(alpha,xm,n,drawInterval)

Example of the inverse transformation method (Pareto function distribution)
- alpha shape parameter of the Pareto distribution
- xm minimum x value for the Pareto distribution
- n number of samples
- drawInterval number of samples after which the histogramms are drawn
"""
function samplePareto(alpha = 1, xm = 1, n = 1000, drawInterval = 50)

    binsuniform = range(0, 1, length = 20)
    binspareto = range(xm-1, xm + 15, length = 100)
    uniformSamples = zeros(n)
    paretoSamples = zeros(n)
    d_uniform = Uniform(0, 1)
    d_pareto = Pareto(alpha, xm)
    anim = Animation()
    @progress for i = 1:n
        newUniformSample = rand(d_uniform)
        # The inverse of the cdf:
        newParetoSample = inverseparetocdf(newUniformSample, alpha, xm)
        uniformSamples[i] = newUniformSample
        paretoSamples[i] = newParetoSample
        if i % min(drawInterval,10^(floor(log10(i)))) == 0
            p1 = histogram(uniformSamples[1:i], bins = binsuniform, normalize = true, label="", lw =0.0, title = "N = $i")
            plot!(x->pdf(d_uniform, x), label= "U(0,1)", lw = 3.0)
            p2 = histogram(paretoSamples[1:i], bins = binspareto, normalize = true, label="", lw = 0.0)
            plot!(x->pdf(d_pareto, x), label= "Pareto", lw = 3.0)
            plot(p1, p2)
            frame(anim)
        end
    end
    return gif(anim, fps = 5)
end


samplePareto()
