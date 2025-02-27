function plot_wigner(ρ_t::Vector{<:Operator}; args...)
    wigner_max_and_min = [
        let W = wigner(ρ, range(-5, stop=5, length=100), range(-5, stop=5, length=100))
            (maximum(W), minimum(W))
        end for ρ in ρ_t
    ]
    wigner_max = maximum([max_t for (max_t, min_t) in wigner_max_and_min])
    wigner_min = minimum([min_t for (max_t, min_t) in wigner_max_and_min])

    anim = @animate for ρ in ρ_t
        plot_wigner(ρ; clims=(wigner_min, wigner_max), args...)
    end

    gif(anim, "wigner_animation.gif", fps=10)
    return anim
end

function plot_wigner(ρ::Operator; c=:viridis, xlabel="x", ylabel="p", title="Wigner Function", kwargs...)
    x = y = range(-5, stop=5, length=100)
    W = wigner(ρ, x, y)

    heatmap(x, y, W; c=c, xlabel=xlabel, ylabel=ylabel, title=title, kwargs...)
end