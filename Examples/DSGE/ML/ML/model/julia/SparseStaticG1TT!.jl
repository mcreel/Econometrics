function SparseStaticG1TT!(T::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})
    SparseStaticResidTT!(T, y, x, params)
@inbounds begin
T[5] = get_power_deriv(y[3],params[1],1)
T[6] = get_power_deriv(y[4],1-params[1],1)
end
    return nothing
end

