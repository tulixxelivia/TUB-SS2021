using Distributions, Random
using Plots; pyplot()

function proposalgauss(mean, std, gaussAccepted, gaussSamples, gaussAcceptance, i)
    x = randn() * std + mean
    if x > 0
        gaussAccepted += 1
        gaussSamples[gaussAccepted] = x
    end
    gaussAcceptance[i] = gaussAccepted/i;
    return gaussAccepted
end


function proposalexp(mean, std, expAccepted, expSamples, expAcceptance, i)
    border = - mean / std
    lambda = (border + sqrt(border ^ 2 + 4)) / 2
    x = rand(Exponential(1 / lambda)) + border
    g = exp(-(x-lambda)^2/2)
    if rand() <= g
        expAccepted += 1
        expSamples[expAccepted] = x - border
    end
    expAcceptance[i] = expAccepted/i;
    return expAccepted
end


"""
    truncatedGaussianSampling(mean, std; n)

Sample from a Gaussian which is truncated to give only positive results
 - mean - mean of the Gaussian
 - std - standard deviation of the Gaussian
 - n - number of samples
"""
function truncatedGaussianSampling(mean = -1, std = 1; n = 1000, drawInterval = 50)
    ## Setup
    # number of accepted samples from Gaussian proposal
    gaussAccepted = 0;
    # number of accepted samples from exponential proposal
    expAccepted = 0;
    # acceptance probabilities over n for Gaussian proposal
    gaussAcceptance = zeros(n);
    # acceptance probabilities over n for exponential proposal
    expAcceptance = zeros(n);
    # accepted samples from Gaussian proposal
    gaussSamples = zeros(n);
    # accepted samples from exponential proposal
    expSamples = zeros(n);

    d = truncated(Normal(mean, std), 0.0, Inf)

    bins = range(0, 6*std, length = 100)

    gaussTime = 0;
    expTime = 0;
    anim = Animation()
    # Sampling Loop
    @progress for i = 1:n
        # Gaussian Proposal
        gaussTime += @elapsed gaussAccepted = proposalgauss(mean, std, gaussAccepted, gaussSamples, gaussAcceptance, i)
        expTime += @elapsed expAccepted = proposalexp(mean, std, expAccepted, expSamples, expAcceptance, i)

        # Plotting
        if i % min(drawInterval, 10 ^ (floor(log10(i)))) == 0
            p1 = plot(gaussAcceptance[1:i], title="Gaussian Proposal\n Acceptance", lab = "", ylims = (0, 1))
            p2 = plot(expAcceptance[1:i], title="Exponential Proposal\n Acceptance", lab = "", ylims = (0, 1))
            p3 = histogram(gaussSamples[1:gaussAccepted], normalize = true, bins = bins, title = "Gaussian Proposal\n Sample Histogram", lw = 0.0, lab="")
            plot!(bins, x->pdf(d, x), lab="", lw=3.0)
            p4 = histogram(expSamples[1:expAccepted], normalize = true, bins = bins, title = "Exponential Proposal\n Sample Histogram", lw = 0.0, lab="")
            plot!(bins, x->pdf(d, x), lab="",lw=3.0)
            plot(p1, p2, p3, p4)
            frame(anim)
        end
    end

    println("Time to get $n samples: Gaussian: $(gaussTime/gaussAcceptance[end]) seconds, Exponential $(expTime/expAcceptance[end])")
    return gif(anim, fps = 8)
end

truncatedGaussianSampling(n = 10000)
