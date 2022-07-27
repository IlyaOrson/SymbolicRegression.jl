using IterTools: ncycle

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
    return [stoic * r for stoic in scoeff]
end

# datasets = Array{Float32}(undef, num_states, num_datasets * num_timepoints)
datasets = []
for ini in initial_conditions
    prob = ODEProblem(f, ini, tspan)
    sol = solve(prob, Rodas4(); saveat=times_per_dataset)
    push!(datasets, Arrat(sol))
end

X = hcat(datasets...)
times = ncycle(times_per_dataset, num_datasets) |> collect
experiments = vcat([fill(Float32(i), num_timepoints) for i in 1:num_datasets]...)
