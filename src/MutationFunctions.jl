module MutationFunctionsModule

import ..CoreModule: CONST_TYPE, Node, copyNode, Options
import ..EquationUtilsModule: countNodes, countConstants, countOperators, countDepth

# Return a random node from the tree
function randomNode(tree::Node)::Node
    if tree.degree == 0
        return tree
    end
    a = countNodes(tree)
    b = 0
    c = 0
    if tree.degree >= 1
        b = countNodes(tree.l)
    end
    if tree.degree == 2
        c = countNodes(tree.r)
    end

    i = rand(1:1+b+c)
    if i <= b
        return randomNode(tree.l)
    elseif i == b + 1
        return tree
    end

    return randomNode(tree.r)
end

# Randomly convert an operator into another one (binary->binary;
# unary->unary)
function mutateOperator(tree::Node, options::Options)::Node
    if countOperators(tree) == 0
        return tree
    end
    node = randomNode(tree)
    while node.degree == 0
        node = randomNode(tree)
    end
    if node.degree == 1
        node.op = rand(1:options.nuna)
    else
        node.op = rand(1:options.nbin)
    end
    return tree
end

# Randomly perturb a constant
function mutateConstant(
        tree::Node, temperature::T,
        options::Options)::Node where {T<:Real}
    # T is between 0 and 1.

    if countConstants(tree) == 0
        return tree
    end
    node = randomNode(tree)
    while node.degree != 0 || node.constant == false
        node = randomNode(tree)
    end

    bottom = convert(T, 1//10)
    maxChange = options.perturbationFactor * temperature + convert(T, 1) + bottom
    factor = maxChange^Float32(rand())
    makeConstBigger = rand() > 0.5

    if makeConstBigger
        node.val *= convert(CONST_TYPE, factor)
    else
        node.val /= convert(CONST_TYPE, factor)
    end

    if rand() > options.probNegate
        node.val *= -1
    end

    return tree
end

# Add a random unary/binary operation to the end of a tree
function appendRandomOp(tree::Node, options::Options, nfeatures::Int; makeNewBinOp::Union{Bool,Nothing}=nothing)::Node
    node = randomNode(tree)
    while node.degree != 0
        node = randomNode(tree)
    end


    if makeNewBinOp === nothing
        choice = rand()
        makeNewBinOp = choice < options.nbin/(options.nuna + options.nbin)
    end

    if makeNewBinOp
        newnode = Node(
            rand(1:options.nbin),
            makeRandomLeaf(nfeatures),
            makeRandomLeaf(nfeatures)
        )
    else
        newnode = Node(
            rand(1:options.nuna),
            makeRandomLeaf(nfeatures)
        )
    end

    if newnode.degree == 2
        node.r = newnode.r
    end
    node.l = newnode.l
    node.op = newnode.op
    node.degree = newnode.degree
    node.val = newnode.val
    node.feature = newnode.feature
    node.constant = newnode.constant

    return tree
end

# Insert random node
function insertRandomOp(tree::Node, options::Options, nfeatures::Int)::Node
    node = randomNode(tree)
    choice = rand()
    makeNewBinOp = choice < options.nbin/(options.nuna + options.nbin)
    left = copyNode(node)

    if makeNewBinOp
        right = makeRandomLeaf(nfeatures)
        newnode = Node(
            rand(1:options.nbin),
            left,
            right
        )
    else
        newnode = Node(
            rand(1:options.nuna),
            left
        )
    end
    if newnode.degree == 2
        node.r = newnode.r
    end
    node.l = newnode.l
    node.op = newnode.op
    node.degree = newnode.degree
    node.val = newnode.val
    node.feature = newnode.feature
    node.constant = newnode.constant
    return tree
end

# Add random node to the top of a tree
function prependRandomOp(tree::Node, options::Options, nfeatures::Int)::Node
    node = tree
    choice = rand()
    makeNewBinOp = choice < options.nbin/(options.nuna + options.nbin)
    left = copyNode(tree)

    if makeNewBinOp
        right = makeRandomLeaf(nfeatures)
        newnode = Node(
            rand(1:options.nbin),
            left,
            right
        )
    else
        newnode = Node(
            rand(1:options.nuna),
            left
        )
    end
    if newnode.degree == 2
        node.r = newnode.r
    end
    node.l = newnode.l
    node.op = newnode.op
    node.degree = newnode.degree
    node.val = newnode.val
    node.feature = newnode.feature
    node.constant = newnode.constant
    return node
end

function makeRandomLeaf(nfeatures::Int)::Node
    if rand() > 0.5
        return Node(randn(CONST_TYPE))
    else
        return Node(rand(1:nfeatures))
    end
end


# Return a random node from the tree with parent, and side ('n' for no parent)
function randomNodeAndParent(tree::Node, parent::Union{Node, Nothing}; side::Char)::Tuple{Node, Union{Node, Nothing}, Char}
    if tree.degree == 0
        return tree, parent, side
    end
    a = countNodes(tree)
    b = 0
    c = 0
    if tree.degree >= 1
        b = countNodes(tree.l)
    end
    if tree.degree == 2
        c = countNodes(tree.r)
    end

    i = rand(1:1+b+c)
    if i <= b
        return randomNodeAndParent(tree.l, tree; side='l')
    elseif i == b + 1
        return tree, parent, side
    end

    return randomNodeAndParent(tree.r, tree; side='r')
end

function randomNodeAndParent(tree::Node)::Tuple{Node, Union{Node, Nothing}, Char}
    return randomNodeAndParent(tree, nothing; side='n')
end

# Select a random node, and replace it an the subtree
# with a variable or constant
function deleteRandomOp(tree::Node, options::Options, nfeatures::Int)::Node
    node, parent, side = randomNodeAndParent(tree)
    isroot = (parent === nothing)

    if node.degree == 0
        # Replace with new constant
        newnode = makeRandomLeaf(nfeatures)
        node.degree = newnode.degree
        node.val = newnode.val
        node.constant = newnode.constant
        if !newnode.constant
            node.feature = newnode.feature
        end
    elseif node.degree == 1
        # Join one of the children with the parent
        if isroot
            return node.l
        elseif parent.l == node
            parent.l = node.l
        else
            parent.r = node.l
        end
    else
        # Join one of the children with the parent
        if rand() < 0.5
            if isroot
                return node.l
            elseif parent.l == node
                parent.l = node.l
            else
                parent.r = node.l
            end
        else
            if isroot
                return node.r
            elseif parent.l == node
                parent.l = node.r
            else
                parent.r = node.r
            end
        end
    end
    return tree
end

# Create a random equation by appending random operators
function genRandomTree(length::Int, options::Options, nfeatures::Int)::Node
    # Note that this base tree is just a placeholder; it will be replaced.
    tree = Node(convert(CONST_TYPE, 1))
    for i=1:length
        # TODO: This can be larger number of nodes than length.
        tree = appendRandomOp(tree, options, nfeatures)
    end
    return tree
end

function genRandomTreeFixedSize(node_count::Int, options::Options, nfeatures::Int)::Node
    tree = makeRandomLeaf(nfeatures)
    cur_size = countNodes(tree)
    while cur_size < node_count
        if cur_size == node_count - 1  # only unary operator allowed.
            options.nuna == 0 && break # We will go over the requested amount, so we must break.
            tree = appendRandomOp(tree, options, nfeatures; makeNewBinOp=false)
        else
            tree = appendRandomOp(tree, options, nfeatures)
        end
        cur_size = countNodes(tree)
    end
    return tree
end

"""Crossover between two expressions"""
function crossoverTrees(tree1::Node, tree2::Node)::Tuple{Node, Node}
    tree1 = copyNode(tree1)
    tree2 = copyNode(tree2)

    node1, parent1, side1 = randomNodeAndParent(tree1)
    node2, parent2, side2 = randomNodeAndParent(tree2)

    node1 = copyNode(node1)

    if side1 == 'l'
        parent1.l = copyNode(node2)
        # tree1 now contains this.
    elseif side1 == 'r'
        parent1.r = copyNode(node2)
        # tree1 now contains this.
    else # 'n'
        # This means that there is no parent2.
        tree1 = copyNode(node2)
    end

    if side2 == 'l'
        parent2.l = node1
    elseif side2 == 'r'
        parent2.r = node1
    else # 'n'
        tree2 = node1
    end
    return tree1, tree2
end

end
