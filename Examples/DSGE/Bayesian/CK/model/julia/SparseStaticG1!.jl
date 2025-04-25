function SparseStaticG1!(T::Vector{<: Real}, g1_v::Vector{<: Real}, y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})
    @assert length(T) >= 6
    @assert length(g1_v) == 28
    @assert length(y) == 11
    @assert length(x) == 2
    @assert length(params) == 15
@inbounds begin
g1_v[1]=1;
g1_v[2]=(-1);
g1_v[3]=(-(get_power_deriv(y[2],(-params[4]),1)));
g1_v[4]=1;
g1_v[5]=(-(T[2]*params[1]*exp(y[6])*get_power_deriv(y[3],params[1]-1,1)));
g1_v[6]=(-(T[4]*exp(y[6])*(1-params[1])*T[5]));
g1_v[7]=(-(T[2]*exp(y[6])*T[5]));
g1_v[8]=1-(1-params[3]);
g1_v[9]=(-(T[1]*T[6]));
g1_v[10]=(-(exp(y[6])*(1-params[1])*T[3]*get_power_deriv(y[4],(-params[1]),1)));
g1_v[11]=(-(exp(y[6])*T[3]*T[6]));
g1_v[12]=1;
g1_v[13]=(-1);
g1_v[14]=(-(T[1]*T[2]));
g1_v[15]=(-(exp(y[6])*(1-params[1])*T[3]*T[4]));
g1_v[16]=1-params[6];
g1_v[17]=(-(T[2]*exp(y[6])*T[3]));
g1_v[18]=(-(params[10]*exp(y[7])));
g1_v[19]=1-params[8];
g1_v[20]=1;
g1_v[21]=1-params[2]*(1+y[10]-params[3]);
g1_v[22]=(-y[9])/(y[8]*y[8]);
g1_v[23]=1;
g1_v[24]=1/y[8];
g1_v[25]=1;
g1_v[26]=(-(y[8]*params[2]));
g1_v[27]=1;
g1_v[28]=(-1);
end
    if ~isreal(g1_v)
        g1_v = real(g1_v)+2*imag(g1_v);
    end
    return nothing
end

