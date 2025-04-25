function SparseDynamicG1!(T::Vector{<: Real}, g1_v::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, steady_state::Vector{<: Real})
    @assert length(T) >= 6
    @assert length(g1_v) == 34
    @assert length(y) == 33
    @assert length(x) == 2
    @assert length(params) == 15
@inbounds begin
g1_v[1]=(-(1-params[3]));
g1_v[2]=(-1);
g1_v[3]=(-params[6]);
g1_v[4]=(-params[8]);
g1_v[5]=1;
g1_v[6]=(-1);
g1_v[7]=(-(get_power_deriv(y[13],(-params[4]),1)));
g1_v[8]=1;
g1_v[9]=(-(T[2]*params[1]*exp(y[17])*get_power_deriv(y[14],params[1]-1,1)));
g1_v[10]=(-(T[4]*exp(y[17])*(1-params[1])*T[5]));
g1_v[11]=(-(T[2]*exp(y[17])*T[5]));
g1_v[12]=1;
g1_v[13]=(-(T[1]*T[6]));
g1_v[14]=(-(exp(y[17])*(1-params[1])*T[3]*get_power_deriv(y[15],(-params[1]),1)));
g1_v[15]=(-(exp(y[17])*T[3]*T[6]));
g1_v[16]=1;
g1_v[17]=(-(T[1]*T[2]));
g1_v[18]=(-(exp(y[17])*(1-params[1])*T[3]*T[4]));
g1_v[19]=1;
g1_v[20]=(-(T[2]*exp(y[17])*T[3]));
g1_v[21]=(-(params[10]*exp(y[18])));
g1_v[22]=1;
g1_v[23]=1;
g1_v[24]=1;
g1_v[25]=(-y[20])/(y[19]*y[19]);
g1_v[26]=1;
g1_v[27]=1/y[19];
g1_v[28]=1;
g1_v[29]=1;
g1_v[30]=(-1);
g1_v[31]=(-(params[2]*(1+y[32]-params[3])));
g1_v[32]=(-(params[2]*y[30]));
g1_v[33]=(-params[7]);
g1_v[34]=(-params[9]);
end
    return nothing
end

