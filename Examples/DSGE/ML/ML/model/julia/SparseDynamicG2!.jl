function SparseDynamicG2!(T::Vector{<: Real}, g2_v::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    @assert length(T) >= 10
    @assert length(g2_v) == 23
    @assert length(y) == 33
    @assert length(x) == 2
    @assert length(params) == 15
@inbounds begin
g2_v[1]=(-(get_power_deriv(y[13],(-params[4]),2)));
g2_v[2]=(-(params[10]*exp(y[18])));
g2_v[3]=(-(T[2]*params[1]*exp(y[17])*get_power_deriv(y[14],params[1]-1,2)));
g2_v[4]=(-(T[5]*T[7]));
g2_v[5]=(-(T[2]*T[5]));
g2_v[6]=(-(T[1]*T[9]));
g2_v[7]=(-(T[1]*T[7]));
g2_v[8]=(-(T[1]*T[2]));
g2_v[9]=(-(T[4]*exp(y[17])*(1-params[1])*T[10]));
g2_v[10]=(-(exp(y[17])*(1-params[1])*T[6]*T[8]));
g2_v[11]=(-(T[4]*exp(y[17])*(1-params[1])*T[6]));
g2_v[12]=(-(exp(y[17])*(1-params[1])*T[3]*get_power_deriv(y[15],(-params[1]),2)));
g2_v[13]=(-(exp(y[17])*(1-params[1])*T[3]*T[8]));
g2_v[14]=(-(exp(y[17])*(1-params[1])*T[3]*T[4]));
g2_v[15]=(-params[2]);
g2_v[16]=(-((-y[20])*(y[19]+y[19])))/(y[19]*y[19]*y[19]*y[19]);
g2_v[17]=(-1)/(y[19]*y[19]);
g2_v[18]=(-(T[2]*exp(y[17])*T[10]));
g2_v[19]=(-(exp(y[17])*T[6]*T[7]));
g2_v[20]=(-(T[2]*exp(y[17])*T[6]));
g2_v[21]=(-(exp(y[17])*T[3]*T[9]));
g2_v[22]=(-(exp(y[17])*T[3]*T[7]));
g2_v[23]=(-(T[2]*exp(y[17])*T[3]));
end
    return nothing
end

