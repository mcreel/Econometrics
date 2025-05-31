"""
    tabular2conv(X::AbstractArray{T, 3}) where T

Transform a (`K × S × T`) array to a CNN format (`1 × T × K × S` array).
`K` is the number of features, `S` is the number of samples, and `T` is the number
of time steps / observations in each sample.

# Arguments
- `X::AbstractArray{T, 3}`: Array to transform.

# Returns
- `Array{T, 4}`: CNN format of `X`.
"""
@views tabular2conv(X) = permutedims(reshape(X, size(X)..., 1), (4, 3, 1, 2))

# In the following losses, Ŷ is the same dimension as Y (CNN style prediction)
rmse_conv(Ŷ, Y) = mean(sqrt.(mean(abs2.(Ŷ - Y), dims=2)))
mse_conv(Ŷ, Y) = mean(abs2, Ŷ - Y)
#export build_tcn, receptive_field_size, necessary_layers

"""
    receptive_field_size(dilation, kernel_size, layers)

Computes the receptive field size for a specified dilation, kernel size, and 
number of layers.

# Arguments
- `dilation::Int`: Dilation factor.
- `kernel_size::Int`: Size of the convolutional kernel.
- `layers::Int`: Number of layers.

# Returns
- `Int`: Receptive field size.
"""
receptive_field_size(dilation::Int, kernel_size::Int, layers::Int) = 
    1 + (kernel_size - 1) * (dilation ^ layers - 1) / (dilation - 1)

"""
    necessary_layers(dilation, kernel_size, receptive_field)

Computes the minimum number of layers necessary to achieve a specified receptive
field size.

# Arguments
- `dilation::Int`: Dilation factor.
- `kernel_size::Int`: Size of the convolutional kernel.
- `receptive_field::Int`: Desired receptive field size.

# Keyword Arguments
- `use_ceil::Bool`: Whether to use `ceil(Int, x)` for the output. (default: `false`)
- `use_floor::Bool`: Whether to use `floor(Int, x)` for the output. (default: `false`)

# Returns
- Minimum number of layers (without `use_ceil` or `use_floor`, this is a `Float64`).
"""
function necessary_layers(
    dilation::Int, kernel_size::Int, receptive_field::Int;
    use_ceil::Bool=false, use_floor::Bool=false
)
    use_ceil && use_floor && throw(
        ArgumentError("Cannot use both `use_ceil` and `use_floor`"))
    
    n = log(dilation, (receptive_field - 1) * (dilation - 1) / (kernel_size - 1)) + 1
    use_ceil ? ceil(Int, n) : use_floor ? floor(Int, n) : n
end

"""
    build_tcn(d::AbstractDGP; dilation_factor, kernel_size, channels, 
        summary_size, dropout_rate, residual, pad)
    
Builds a TCN for a specified DGP.

# Arguments
- `d::AbstractDGP`: DGP to build the TCN for.

# Keyword Arguments
- `nlayers::Int`: Number of layers. (default: `0`, computes the necessary 
    layers using [`receptive_field_size`](@ref))
- `dilation_factor::Int`: Dilation factor. (default: `2`)
- `kernel_size::Int`: Size of the convolutional kernel. (default: `8`)
- `channels::Int`: Number of channels. (default: `32`)
- `receptive_field_size::Int`: Receptive field size. (default: `0`, uses the 
    length of the [`AbstractDGP`](@ref) instead)
- `dropout_rate::AbstractFloat`: Dropout rate. (default: `0.2`)
- `summary_size::Int`: Size of the summary convolution. (default: `10`)
- `residual::Bool`: Whether to use residual connections. (default: `true`)
- `pad`: Padding function. (default: `SamePad()`)
- `dev::Function`: Device to use. (default: `cpu`)

# Returns
- `Chain`: TCN for the specified DGP.
"""
function build_tcn(N, dim_in, dim_out; 
    nlayers::Int=0, dilation_factor::Int=2, kernel_size::Int=8, channels::Int=32, 
    receptive_field_size::Int=160, dropout_rate::AbstractFloat=0., summary_size::Int=10, 
    residual::Bool=true,  pad=SamePad(), dev::Function=cpu
)

    if nlayers == 0
        nlayers = necessary_layers(dilation_factor, kernel_size, receptive_field_size,
            use_ceil=true)
    end

    dev(
        Chain(
            TCN(
                vcat(dim_in, [channels for _ in 1:nlayers], 1),
                kernel_size=kernel_size, dropout_rate=dropout_rate, pad=pad,
                residual=residual, dilation_factor=dilation_factor
            ),
            Conv((1, summary_size), 1 => 1, stride=summary_size),
            Flux.flatten,
            Dense(N ÷ summary_size => N ÷ summary_size, hardtanh), # this is a new layer
            Dense(N ÷ summary_size => dim_out)
        )
    )
end

"""
    TemporalBlock(chan_in, chan_out; dilation, kernel_size, residual, pad, dropout_rate)

Temporal block with `chan_in` input channels and `chan_out` output channels. Each
block consists of two causal convolutional layers with `kernel_size` and `dilation`
followed by batch normalization and dropout. If `residual` is `true`, a skip
connection is added.

# Arguments
- `chan_in::Int`: Number of input channels.
- `chan_out::Int`: Number of output channels.

# Keyword arguments
- `dilation::Int`: Kernel dilation.
- `kernel_size::Int`: Size of the convolutional kernel.
- `residual::Bool`: Whether to use residual connections.
- `pad`: Padding to use for the convolutional layers.
- `dropout_rate::AbstractFloat`: Dropout rate to use for the convolutional layers.

# Returns
- `Chain`: Temporal block
"""
function TemporalBlock(
    chan_in::Int, chan_out::Int; 
    dilation::Int, kernel_size::Int,
    residual::Bool, pad, dropout_rate::AbstractFloat
)
    # Causal convolutions
    causal_conv = Chain(
        Conv((1, kernel_size), chan_in => chan_out, dilation = dilation, 
            pad = pad),
        BatchNorm(chan_out, leakyrelu),
        Dropout(dropout_rate),
        Conv((1, kernel_size), chan_out => chan_out, dilation = dilation, 
            pad = pad),
        BatchNorm(chan_out, leakyrelu),
        Dropout(dropout_rate),
    )
    residual || return causal_conv
    # Skip connection (residual net)
    residual_conv = Conv((1, 1), chan_in => chan_out)
    Chain(
        Parallel(+, causal_conv, residual_conv),
        x -> leakyrelu.(x)
    )
end


"""
    TCN(channels; kernel_size, dilation_factor, residual, pad, dropout_rate)

Temporal convolutional network (TCN) with `length(channels) - 1` layers. Each layer
is a `TemporalBlock` with `channels[i]` input channels and `channels[i+1]`

# Arguments
- `channels::AbstractVector{Int}`: Number of input and output channels for each layer.
- `kernel_size::Int`: Size of the convolutional kernel.
- `dilation_factor::Int`: Factor by which the dilation is increased for each layer. (default: `2`)
- `residual::Bool`: Whether to use residual connections. (default: `true`)
- `pad`: Padding to use for the convolutional layers. (default: `SamePad()`)
- `dropout_rate::AbstractFloat`: Dropout rate to use for the convolutional layers. (default: `0.`)

# Returns
- `Chain`: TCN
"""
function TCN(
    channels::AbstractVector{Int}; 
    kernel_size::Int, dilation_factor::Int=2, residual::Bool=true, 
    pad=SamePad(), dropout_rate::AbstractFloat=0.,
)
    Chain([
        TemporalBlock(chan_in, chan_out, dilation=dilation_factor ^ (i - 1), 
            kernel_size=kernel_size, residual=residual, pad=pad,
            dropout_rate=dropout_rate
        ) 
        for (i, (chan_in, chan_out)) ∈ enumerate(zip(channels[1:end-1], channels[2:end]))]...)
end


"""
    train_network!(net::MomentNetwork, dgp::AbstractDGP; verbose::Bool=true)

Train a moment network on a DGP.

# Arguments
- `net::MomentNetwork`: The moment network.
- `dgp::AbstractDGP`: The DGP to train on.
- `verbose::Bool`: Whether to print information about the training process.

# Returns
- `iterations`: The iterations at which the losses were computed.
- `losses`: The losses at each iteration.
"""
function train_network!(net; 
        batchsize=128, 
        validation_size=1000, 
        nsamples=1000, 
        epochs=2,
        verbose::Bool=true,
        validation_freq = 10, # validate every X samples
        print_every = 10      # print to screen every X iterations
)
    Flux.trainmode!(net) # Ensure that the model is in training mode
    θ = Flux.params(net) # Extract model parameters
    
    best_model = deepcopy(net) # Reset the best model
    best_loss = Inf # Initialize the best loss

    validation = validation_size > 0
    loss_iteration = 1

    # Handle validation
    if validation
        Xval, Yval = MakeData(validation_size)
        verbose && @info "Validation set size: $(size(Yval, 2))"
        Xval = gpu(Xval)
        Yval = gpu(Yval)
    end

    if verbose # Compute pre-training loss
        Ŷ = net(Xval)
        loss = rmse_conv(Ŷ, Yval)
        @info "Pre-training loss: $(loss)"
    end

    for sample ∈ 1:nsamples
        # Generate training data
        X, Y = MakeData(batchsize)
        X = gpu(X)
        Y = gpu(Y)

        # Gradient step: repeat with same data epochs times
        for epoch = 1:epochs
            ∇ = gradient(θ) do 
                Ŷ = net(X)
                rmse_conv(Ŷ, Y)
            end
            Flux.update!(ADAMW(0.001), θ, ∇)
        end

        # Handle validation
        if validation && mod(sample, validation_freq) == 0
            Ŷ = net(Xval)
            loss = rmse_conv(Ŷ, Yval)
            if loss < best_loss
                best_loss = loss
                best_model = deepcopy(net)
                printstyled("best so far: $(best_loss)\n", color = :green)

                if sample > 1000
                    net_state = cpu(Flux.state(net))
                    BSON.@save "best_model.bson"  net_state
                end    
            end
            if verbose && mod(sample, print_every) == 0
                @info "Sample: $(sample), loss: $(loss)"
            end
        else
            if verbose && mod(sample, print_every) == 0
                @info "Sample: $(sample)"
            end
        end
    end
    net = deepcopy(best_model)
    nothing
end

