import Pkg
# project_dir = "/Users/md1621/Desktop/PhD-Code/Physics-Informed_ADoK/physics_informed_SR"
# exp_dir = "/Users/md1621/Desktop/PhD-Code/Physics-Informed_ADoK/physics_informed_SR/exp_data"
# rate_dir = "/Users/md1621/Desktop/PhD-Code/Physics-Informed_ADoK/physics_informed_SR/const_data"
# Pkg.activate(project_dir)
# Pkg.instantiate()

# exp_dir = pwd() * "exp_data"
exp_dir = "exp_data"

using IterTools: ncycle
using SymbolicRegression
using Infiltrator
using DelimitedFiles

num_datasets = 5
num_states = 4

tspan = (0e0, 1e1)
num_timepoints = 15

times_per_dataset=collect(range(tspan[begin], tspan[end]; length=num_timepoints))

ini_T = [1e0, 5e0, 5e0, 1e0, 1e0]
ini_H = [8e0, 8e0, 3e0, 3e0, 8e0]
ini_B = [2e0, 0e0, 0e0, 0e0, 2e0]
ini_M = [3e0, 5e-1, 5e-1, 3e0, 5e-1]

# READ THE DATASET FROM PYTHON
for i in 1:num_datasets
    datasets = readdlm(exp_dir*"/exp_$i.csv", ',', Float64, '\n')
    #------------------------------#

    for j in 1:num_states
        X = reshape(times_per_dataset, 1, :)
        y = reshape(datasets[j, :], 1, :)

        if j == 1
            name = "hof_files/hall_of_fame_T$i.csv"
        elseif j == 2
            name = "hof_files/hall_of_fame_H$i.csv"
        elseif j == 3
            name = "hof_files/hall_of_fame_B$i.csv"
        elseif j == 4
            name = "hof_files/hall_of_fame_M$i.csv"
        end

        options = Options(; # NOTE add new constraint here
            binary_operators=[+, *, /, -], unary_operators=[exp], constraint_initial_condition=true, hofFile=name
        )

        hall_of_fame = equation_search(
            X, y, niterations=2, options=options
        )
    end
end

# NOTE model selection done in python to generate the following file...

#=
conc_data = readdlm(rate_dir*"/conc_data_for_rate_models.csv", ',', Float64, '\n')

for j in 1:num_states
    X = reshape(conc_data, num_states, :)

    if j == 1
        name = "hof_files/hall_of_fame_rate_T$i.csv"
        a = readdlm(rate_dir*"/rate_data_T.csv", ',', Float64, '\n')
        y = reshape(a, 1, :)
    elseif j == 2
        name = "hof_files/hall_of_fame_rate_H$i.csv"
        a = readdlm(rate_dir*"/rate_data_H.csv", ',', Float64, '\n')
        y = reshape(a, 1, :)
    elseif j == 3
        name = "hof_files/hall_of_fame_rate_B$i.csv"
        a = readdlm(rate_dir*"/rate_data_B.csv", ',', Float64, '\n')
        y = reshape(a, 1, :)
    elseif j == 4
        name = "hof_files/hall_of_fame_rate_M$i.csv"
        a = readdlm(rate_dir*"/rate_data_M.csv", ',', Float64, '\n')
        y = reshape(a, 1, :)
    end

    options = Options(; # NOTE add new constraint here
        binary_operators=[+, *, /, -], constraint_initial_condition=false, hofFile=name
    )

    hall_of_fame = equation_search(
        X, y, niterations=2, options=options
    )
end
=#

#TODO: make sure the directories are correct (pretty sure they are not)
#TODO: declare the name of the variables for the models for easier manipulation
