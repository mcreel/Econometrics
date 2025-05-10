function SparseStaticResidTT!(T::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})
@inbounds begin
T[1] = params[1]*exp(y[6])*y[3]^(params[1]-1)
T[2] = y[4]^(1-params[1])
T[3] = y[3]^params[1]
T[4] = y[4]^(-params[1])
end
    return nothing
end

