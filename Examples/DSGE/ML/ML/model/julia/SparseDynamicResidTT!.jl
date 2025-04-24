function SparseDynamicResidTT!(T::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
@inbounds begin
T[1] = params[1]*exp(y[17])*y[14]^(params[1]-1)
T[2] = y[15]^(1-params[1])
T[3] = y[14]^params[1]
T[4] = y[15]^(-params[1])
end
    return nothing
end

