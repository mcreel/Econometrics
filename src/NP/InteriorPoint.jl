#  taken from QuantileRegression.jl, then modified a bit
using LinearAlgebra

macro atleastversion(x)
    if eval(x) >= v"0.4-"
        return :(qrfact(Q*X, Val{true}))
    else
        return :(qrfact(Q*X, pivot = true))
    end
end

# A Frisch-Newton algorithm as described in http://projecteuclid.org/euclid.ss/1030037960
# qreg_ip_coef is a modified lp_fnm of Daniel Morillo & Roger Koenker
# Translated from Ox to Matlab by Paul Eilers 1999
# Modified by Roger Koenker 2000
# More changes by Paul Eilers 2004
# QuantileRegression.jl developers 2015
#
# The functions in this file is released under the BSD-3 license, although the original
# versions were not. Permissing was given to Patrick Kofod Mogensen by Roger Koenker and
# Paul Eilers on April 24 2015 via E-mail. Original mail correspondences can be provided
# upon request.

function bound(x, dx)
# Fill vector with allowed step lengths
# Replace with -x/dx for negative dx
    b = 1e20 .+ 0.0 .* x
    for i = 1:length(dx)
        if dx[i] < 0.0
            @inbounds b[i] = -x[i] / dx[i]
        end
    end
    return b
end

function qreg_coef(Y::Vector, X::Matrix, p)
# input: X is an n x k matrix of exogenous regressors,
#        Y is an n x 1 vector of outcome variables
#        p \in (0,1) is the quantile of interest
# Output: p^{th} regression quantiles.
# Construct the dual problem of quantile regression
c = -Y

m = size(X, 1)
u = ones(m)
x = (1 - p) .* u
b = X'x

# Solve a linear program by the interior point method:
# min(c * u), s.t. A * x = b and 0 < x < u
# An initial feasible solution has to be provided as x

# Set some constants
beta = 0.99995
small = 1e-6
max_it = 50
n, m = size(X)

# Generate inital feasible point
s = u - x
y = -X\Y
r = c - X*y
BLAS.axpy!(0.001, (r .== 0.0).*1.0, r)
z = r .* (r .> 0.0)
w = z - r
gap = dot(c, x) - dot(y, b) + dot(w, u)

# Start iterations
it = 0
for it = 1:max_it
    #   Compute affine step
    q = 1 ./ (z ./ x + w ./ s)
    r = z - w
    Q = Diagonal(sqrt.(q)) # Very efficient to do since Q diagonal
    # AQtF = @atleastversion(VERSION)
    AQtF = qr(Q*X) # PE 2004
    rhs = Q*r        # "
    dy = AQtF\rhs   # "
    dx = q.*(X*dy - r)
    ds = -dx
    dz = -z .* (1.0 .+ dx ./ x)
    dw = -w .* (1.0 .+ ds ./ s)

    # Compute maximum allowable step lengths
    fx = bound(x, dx)
    fs = bound(s, ds)
    fw = bound(w, dw)
    fz = bound(z, dz)
    fpv = min.(fx, fs)
    fdv = min.(fw, fz)
    fp = min.(minimum(beta * fpv), 1)
    fd = min.(minimum(beta * fdv), 1)

    # If full step is feasible, take it. Otherwise modify it
    if min.(fp, fd) < 1.0

        # Update mu
        mu = dot(z, x) + dot(w, s)
        g = dot(z + fd*dz, x + fp*dx) + dot(w + fd*dw, s + fp*ds)
        mu = mu * (g / mu)^3 / (2 * n)

        # Compute modified step
        dxdz = dx .* dz
        dsdw = ds .* dw
        xinv = 1 ./ x
        sinv = 1 ./ s
        xi = mu .* (xinv - sinv)
        #rhs = rhs + Q * (dxdz - dsdw - xi)
        BLAS.axpy!(1.0, Q * (dxdz - dsdw - xi), rhs) # no gemv-wrapper gemv(Q, (dxdz - dsdw - xi), rhs,1,1,n)?

        dy = AQtF\rhs
        dx = q .* (X*dy + xi - r - dxdz + dsdw)
        ds = -dx
        for i = 1:length(dz)
            dz[i] = mu * xinv[i] - z[i] - xinv[i] * z[i] * dx[i] - dxdz[i]
            dw[i] = mu * sinv[i] - w[i] - sinv[i] * w[i] * ds[i] - dsdw[i]
        end

        # Compute maximum allowable step lengths
        fx = bound(x, dx)
        fs = bound(s, ds)
        fw = bound(w, dw)
        fz = bound(z, dz)
        fp = min.(fx, fs)
        fd = min.(fw, fz)
        fp = min.(minimum(beta .* fp), 1)
        fd = min.(minimum(beta .* fd), 1)

    end

    # Take the steps
    BLAS.axpy!(fp, dx, x)
    BLAS.axpy!(fp, ds, s)
    BLAS.axpy!(fd, dy, y)
    BLAS.axpy!(fd, dw, w)
    BLAS.axpy!(fd, dz, z)

    gap = dot(c, x) - dot(y, b) + dot(w, u)

    if gap < small
        break
    end
end

return -y
end
