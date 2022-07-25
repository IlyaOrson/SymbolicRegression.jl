module LossFunctionsModule

import Random: randperm
import LossFunctions: value, AggMode, SupervisedLoss
import ..CoreModule: Options, Dataset, Node
import ..EquationUtilsModule: compute_complexity
import ..EvaluateEquationModule: eval_tree_array, differentiable_eval_tree_array

using IterTools: groupby
using Infiltrator: @infiltrate
using OrdinaryDiffEq

function loss( # fmt: off
    x::AbstractArray{T},
    y::AbstractArray{T},
    options::Options{A,B,dA,dB,C,D}, # fmt: on
)::T where {T<:Real,A,B,dA,dB,C,D}
    @error "default loss"
    if C <: SupervisedLoss
        return value(options.loss, y, x, AggMode.Mean())
    elseif C <: Function
        return sum(options.loss.(x, y)) / length(y)
    else
        error("Unrecognized type for loss function: $(C)")
    end
end

function loss(
    x::AbstractArray{T},
    y::AbstractArray{T},
    w::AbstractArray{T},
    options::Options{A,B,dA,dB,C,D},
)::T where {T<:Real,A,B,dA,dB,C,D}
    @error "default loss weighted"
    if C <: SupervisedLoss
        return value(options.loss, y, x, AggMode.WeightedMean(w))
    elseif C <: Function
        return sum(options.loss.(x, y, w)) / sum(w)
    else
        error("Unrecognized type for loss function: $(C)")
    end
end

# Evaluate the loss of a particular expression on the input dataset.
function eval_loss(tree::Node{T}, dataset::Dataset{T}, options::Options)::T where {T<:Real}
    (prediction, completion) = eval_tree_array(tree, dataset.X, options)
    if !completion
        return T(Inf)
    end

    function integrate(initial_condition, tspan)
        function f(u, p, t)
            matrix_state = reshape(collect(u), length(u), :)
            (prediction, completion) = eval_tree_array(tree, matrix_state, options)
            return prediction[begin]
        end
        prob = ODEProblem(f, initial_condition, tspan)
        # local sol
        # try
            sol = solve(prob, Rodas4P(); verbose=false, maxiters=10^3)  # 1.0e-7 default dtmin, 1e5 default maxiters
        # catch
        #     @infiltrate
        # end
        # @infiltrate sol.retcode != :Success
        return sol.u[end], sol.retcode
    end

    # times = dataset.times
    # concentrations = dataset.states

    # timepoints = 21
    # times = range(0f0, 5f0, timepoints)  # FIXME
    # X = Matrix(reshape(times, 1, timepoints))
    # y = exp.(X ./ 2) .+ 5
    # X .= y

    loss_ = 0f0

    # integrate between states
    # for (i, (t, s)) in enumerate(zip(times, states))
    #     i == 1 && continue
    #     timespan = t - times[i-1]
    #     initial_condition = states[i-1]
    #     prediction = integrate(initial_condition, timespan)
    #     loss_ += (prediction - s) .^ 2
    # end

    states = dataset.X
    times = dataset.times
    experiments = dataset.experiments

    @assert length(times) == length(experiments) == size(states, 2)
    per_experiment = groupby(x -> x[begin], zip(experiments, times, eachcol(states)))
    for experiment in per_experiment
        _, t0, initial_condition = experiment[begin]  # supposing it is sorted
        for (i, data) in enumerate(experiment)
            _, t, s = data  # experiment tag is no longer used
            if i == 1
                # integrate from first state always
                continue
            end
            @infiltrate i != 1
            timespan = (t0, t)
            prediction, retcode = integrate(initial_condition, timespan)
            if retcode != :Success
                loss_ = Inf
                break
            end
            loss_ += (prediction - s) .^ 2
        end
    end
    # @info tree, loss_, compute_complexity(tree, options)
    #@infiltrate
    return loss_
    # if dataset.weighted
    #     return loss(prediction, dataset.y, dataset.weights, options)
    # else
    #     return loss(prediction, dataset.y, options)
    # end
end

# Compute a score which includes a complexity penalty in the loss
function loss_to_score(
    loss::T, baseline::T, tree::Node{T}, options::Options
)::T where {T<:Real}
    normalization = if baseline < T(0.01)
        T(0.01)
    else
        baseline
    end
    normalized_loss_term = loss / normalization
    size = compute_complexity(tree, options)
    parsimony_term = size * options.parsimony
    normalized_loss_term + parsimony_term < 0.0
    return normalized_loss_term + parsimony_term
end

# Score an equation
function score_func(
    dataset::Dataset{T}, tree::Node{T}, options::Options
)::Tuple{T,T} where {T<:Real}
    result_loss = eval_loss(tree, dataset, options)
    score = loss_to_score(result_loss, dataset.baseline_loss, tree, options)
    return score, result_loss
end

# Score an equation with a small batch
function score_func_batch(
    dataset::Dataset{T}, tree::Node{T}, options::Options
)::Tuple{T,T} where {T<:Real}
    batch_idx = randperm(dataset.n)[1:(options.batchSize)]
    batch_X = dataset.X[:, batch_idx]
    batch_y = dataset.y[batch_idx]
    (prediction, completion) = eval_tree_array(tree, batch_X, options)
    if !completion
        return T(0), T(Inf)
    end

    if !dataset.weighted
        result_loss = loss(prediction, batch_y, options)
    else
        batch_w = dataset.weights[batch_idx]
        result_loss = loss(prediction, batch_y, batch_w, options)
    end
    score = loss_to_score(result_loss, dataset.baseline_loss, tree, options)
    return score, result_loss
end

function update_baseline_loss!(dataset::Dataset{T}, options::Options) where {T<:Real}
    dataset.baseline_loss = one(T)  # [1f0]
    # dataset.baseline_loss = if dataset.weighted
    #     loss(dataset.y, ones(T, dataset.n) .* dataset.avg_y, dataset.weights, options)
    # else
    #     loss(dataset.y, ones(T, dataset.n) .* dataset.avg_y, options)
    # end
    return nothing
end

end
