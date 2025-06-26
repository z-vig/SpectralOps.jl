# utils/linear_interpolation.jl

"
    linear_interpolation(x_pts, y_pts, interp_x)

Perform a linear interpolation between x_pts and fill 
"
function linear_interpolation(
    x_pts::AbstractVector{T1},
    y_pts::AbstractVector{T2},
    interp_x::AbstractVector{T3}
)::Vector{Float64} where {T1<:Real, T2<:Real, T3<:Real}
    interp = Array{Float64}(undef, length(interp_x))
    slopes = Array{Float64}(undef, length(x_pts) - 1)
    intercepts = Array{Float64}(undef, length(x_pts) - 1)

    # Filling in between points.
    for n ∈ range(1, length(x_pts)-1)
        x1 = x_pts[n]
        x2 = x_pts[n+1]
        y1 = y_pts[n]
        y2 = y_pts[n+1]
        m = (y2 - y1) / (x2 - x1)
        b = y2 - m * x2
        slopes[n] = m
        intercepts[n] = b

        idx = findall(x -> (x <= x2) & (x > x1), interp_x)
        interp[idx] = m .* interp_x[idx] .+ b
    end

    # Extrapolating where interp_x is greater than or less than points.
    for n in findall(x -> x <= x_pts[1], interp_x)
        interp[n] = interp_x[n] * slopes[1] + intercepts[1]
    end

    for n in findall(x -> x > x_pts[end], interp_x)
        interp[n] = interp_x[n] * slopes[end] + intercepts[end]
    end

    return interp

end