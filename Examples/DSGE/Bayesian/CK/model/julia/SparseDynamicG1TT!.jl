function SparseDynamicG1TT!(T::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    SparseDynamicResidTT!(T, y, x, params, steady_state)
@inbounds begin
T[5] = get_power_deriv(y[14],params[1],1)
T[6] = get_power_deriv(y[15],1-params[1],1)
end
    return nothing
end

