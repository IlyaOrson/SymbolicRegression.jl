using IterTools: ncycle
using OrdinaryDiffEq

using SymbolicRegression

num_datasets = 6
num_timepoints = 10
num_states = 2

tspan = (0f0,10f0)
times_per_dataset=Float32.(collect(range(tspan[begin], tspan[end], num_timepoints)))

ini_Ca = range(2f0,10f0, num_datasets)
ini_Cb = range(0f0,2f0, num_datasets)
initial_conditions = [[x0[begin],x0[end]] for x0 in zip(ini_Ca, ini_Cb)]

scoeff = [-1f0, 1f0]

function rate(Ca,Cb; kf=7f0, kr=3f0, ka=4f0, kb=2f0, kc=6f0)
    num = kf*Ca - kr*Cb
    den = ka*Ca + kb*Cb + kc
    return num / den
end

function f(u,p,t)
    Ca, Cb = u
    r = rate(Ca, Cb)
    return [(1 / stoic) * r for stoic in scoeff]
end

function generate_datasets(; noise_per_concentration=nothing)
    datasets = []
    for ini in initial_conditions
        prob = ODEProblem(f, ini, tspan)
        sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=times_per_dataset)
        arr = Array(sol)
        if isnothing(noise_per_concentration)
            push!(datasets, Array(sol))
        else
            noise_matrix = vcat([noise_level * randn(Float32, (1,length(times_per_dataset))) for noise_level in noise_per_concentration]...)
            push!(datasets, Array(sol) .+ noise_matrix)
            # push!(datasets, Array(sol))
        end
    end
    return datasets
end

datasets = generate_datasets(; noise_per_concentration=[0.13552795534109782f0, 0.2144720446589021f0])
X = hcat(datasets...)
times = ncycle(times_per_dataset, num_datasets) |> collect
experiments = vcat([fill(Float32(i), num_timepoints) for i in 1:num_datasets]...)

y = X[1,:]

options = SymbolicRegression.Options(binary_operators=(+, *, /, -))
hall_of_fame = EquationSearch(X, y, niterations=40, options=options, numprocs=4, times=times, experiments=experiments, stoic_coeff=scoeff)

dominating = calculate_pareto_frontier(X, y, hall_of_fame, options)

println("Complexity\tMSE\tEquation")

for member in dominating
    complexity = compute_complexity(member.tree, options)
    loss = member.loss
    string = string_tree(member.tree, options)

    println("$(complexity)\t$(loss)\t$(string)")
end

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
