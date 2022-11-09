import Pkg
# project_dir = "/rds/general/user/os220/home/SymbolicRegression.jl"
# Pkg.activate(project_dir)
# Pkg.instantiate()

using IterTools: ncycle
using OrdinaryDiffEq
using SymbolicRegression
using Infiltrator
using DelimitedFiles

num_datasets = 1
num_states = 3

scoeff = [-1e0, -3e0, 1e0]

# READ THE DATASET FROM PYTHON
# num_initial_conditions = 3
# datasets = [permutedims(readdlm(project_dir*"/data_$i.csv", '|', Float64, '\n')) for i in 0:num_initial_conditions-1]

# X = hcat(datasets...)
# times = ncycle(times_per_dataset, num_datasets) |> collect
# experiments = vcat([fill(Float64(i), num_timepoints) for i in 1:num_datasets]...)

# y = X[1,:]

# NOTE in hpc add project_dir *
X = readdlm("concentration_data.csv", '|', Float64, '\n')
y = readdlm("rate_data.csv", '|', Float64, '\n')
times = nothing
experiments = nothing

options = SymbolicRegression.Options(binary_operators=(+, *, /, -))
hall_of_fame = EquationSearch(X, y, niterations=40, options=options, numprocs=8, times=times, experiments=experiments, stoic_coeff=scoeff)

dominating = calculate_pareto_frontier(X, y, hall_of_fame, options)

println("Complexity\tMSE\tEquation")

for member in dominating
    complexity = compute_complexity(member.tree, options)
    loss = member.loss
    string = string_tree(member.tree, options)

    println("$(complexity)\t$(loss)\t$(string)")
end

#=
# test against symbolic solution
proposed_rate(x1,x2) = ((x1 - ((x2 - x1) / 1.3333641f0)) / ((((x2 - -0.15033427f0) * 1.500032f0) - x2) + (x1 + 1.2743267f0)))
function f_(u,p,t)
    Ca, Cb = u
    r = proposed_rate(Ca, Cb)
    return [stoic * r for stoic in scoeff]
end
datasets_ = []
for ini in initial_conditions
    prob = ODEProblem(f_, ini, tspan)
    sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=times_per_dataset)
    push!(datasets_, Array(sol))
end
=#
