include("test_params.jl")
using SymbolicRegression, Test

n = 10

options = Options(;
    default_params...,
    binary_operators=(+, -, *, /),
    unary_operators=(cos, sin),
    probPickFirst=0.999,
    ns=n,
)

for reverse in [false, true]
    members = PopMember{Float32}[]

    # Generate members with scores from 0 to 1:
    for i=1:n
        tree = Node("x1") * 3.2f0
        score = Float32(i-1)/(n-1)
        if reverse
            score = 1 - score
        end
        loss = 1f0  # (arbitrary for this test)
        push!(members, PopMember(tree, score, loss))
    end

    pop = Population(members, n)

    dummy_frequencies = [0f-10 for i=1:100]
    best_pop_member = [
        SymbolicRegression.bestOfSample(pop, dummy_frequencies, options).score
        for j=1:100
    ]

    mean_value = sum(best_pop_member)/length(best_pop_member)

    # Make sure average score is small
    @test mean_value < 0.1
end