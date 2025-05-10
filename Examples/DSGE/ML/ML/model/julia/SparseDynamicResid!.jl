function SparseDynamicResid!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    @assert length(T) >= 4
    @assert length(residual) == 11
    @assert length(y) == 33
    @assert length(x) == 2
    @assert length(params) == 15
@inbounds begin
    residual[1] = (y[19]) - (y[13]^(-params[4]));
    residual[2] = (y[20]) - (params[10]*exp(y[18]));
    residual[3] = (y[21]) - (T[1]*T[2]);
    residual[4] = (y[22]) - (exp(y[17])*(1-params[1])*T[3]*T[4]);
    residual[5] = (y[19]) - (params[2]*y[30]*(1+y[32]-params[3]));
    residual[6] = (y[20]/y[19]) - (y[22]);
    residual[7] = (y[17]) - (params[6]*y[6]+params[7]*x[1]);
    residual[8] = (y[18]) - (params[8]*y[7]+params[9]*x[2]);
    residual[9] = (y[12]) - (T[2]*exp(y[17])*T[3]);
    residual[10] = (y[16]) - (y[12]-y[13]);
    residual[11] = (y[14]) - (y[5]+(1-params[3])*y[3]);
end
    return nothing
end

