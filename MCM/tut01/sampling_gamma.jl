using Distributions, Random
using SpecialFunctions
using Plots; pyplot()

function proposalcauchy(a, k, accepted, samples, acceptance, i)
    proposal = sqrt(2*a-1)*tan(pi*rand()-0.5*π)+(a-1);
    u = rand()*k/(1+((proposal-(a-1))^2)/(2*a-1));
    if (proposal >= 0) && (u < (proposal ^ (a - 1)) * exp(-proposal) / gamma(a))
        accepted += 1
        samples[accepted] = proposal
    end
    acceptance[i] = accepted/i;
    return accepted
end


"""
    rejectionSampleGamma(a, std; n)

Sample from a gamma distribution by a rejection sampler with a cauchy distribution as proposal
 - a - Shape parameter of the gamma
 - n - number of samples
"""
function rejectionSampleGamma(a = 2; n = 1000, drawInterval = 50)
    ## Setup
    # number of accepted samples
    accepted = 0;
    # acceptance probabilities over n
    acceptance = zeros(n);
    # accepted samples
    samples = zeros(n);

    b = 1.0
    d = Gamma(a, b)
    d_proposal = Cauchy()
    k = ((a-1)^(a-1))*exp(-(a-1))/gamma(a)
    bins = range(0, 10, length = 100)

    time = 0.0;
    anim = Animation()
    # Sampling Loop
    @progress for i = 1:n
        # Gaussian Proposal
        time += @elapsed accepted = proposalcauchy(a, k, accepted, samples, acceptance, i)

        # Plotting
        if i % min(drawInterval, 10 ^ (floor(log10(i)))) == 0
            p1 = plot(acceptance[1:i], title="Proposal Acceptance", lab="", ylims = (0, 1))
            p2 = histogram(samples[1:accepted], normalize = true, bins = bins, title = "Proposal Sample\n Histogram", lw = 0.0, lab="")
            plot!(bins, x->pdf(d, x), lab="", lw=3.0)
            plot(p1, p2)
            frame(anim)
        end
    end

    println("Time to get $n samples: $(time/acceptance[end]) seconds")
    return gif(anim, fps = 8)
end

rejectionSampleGamma(n = 10000)
