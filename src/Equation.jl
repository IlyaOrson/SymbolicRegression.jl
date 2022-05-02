module EquationModule

import ..ProgramConstantsModule: CONST_TYPE
import ..OptionsStructModule: Options

# Define a serialization format for the symbolic equations:
mutable struct Node
    #Holds operators, variables, constants in a tree
    degree::Int #0 for constant/variable, 1 for cos/sin, 2 for +/* etc.
    constant::Bool #false if variable
    val::CONST_TYPE
    # ------------------- (possibly undefined below)
    feature::Int #Either const value, or enumerates variable.
    op::Int #enumerates operator (separately for degree=1,2)
    l::Node
    r::Node

    Node(d::Int, c::Bool, v::CONST_TYPE) = new(d, c, v)
    Node(d::Int, c::Bool, v::CONST_TYPE, f::Int) = new(d, c, v, f)
    Node(d::Int, c::Bool, v::CONST_TYPE, f::Int, o::Int, l::Node) = new(d, c, v, f, o, l)
    Node(d::Int, c::Bool, v::CONST_TYPE, f::Int, o::Int, l::Node, r::Node) = new(d, c, v, f, o, l, r)
end

Node(val::CONST_TYPE) =                                                     Node(0, true,                       val                                     ) #Leave other values undefined
"""
    Node(feature::Int)

Create a variable node using feature `feature::Int`
"""
Node(feature::Int) =                                                        Node(0, false, convert(CONST_TYPE, 0f0), feature                            )
"""
    Node(op::Int, l::Node)

Apply unary operator `op` (enumerating over the order given) to `Node` `l`
"""
Node(op::Int, l::Node) =                                                    Node(1, false, convert(CONST_TYPE, 0f0),       0,      op,        l         )
"""
    Node(op::Int, l::Union{AbstractFloat, Int})

Short-form for creating a scalar/variable node, and applying a unary operator
"""
Node(op::Int, l::Union{AbstractFloat, Int}) =                               Node(1, false, convert(CONST_TYPE, 0f0),       0,      op,  Node(l)         )
"""
    Node(op::Int, l::Node, r::Node)

Apply binary operator `op` (enumerating over the order given) to `Node`s `l` and `r`
"""
Node(op::Int, l::Node, r::Node) =                                           Node(2, false, convert(CONST_TYPE, 0f0),       0,      op,        l,       r)
"""
    Node(op::Int, l::Union{AbstractFloat, Int}, r::Node)

Short-form to create a scalar/variable node, and apply a binary operator
"""
Node(op::Int, l::Union{AbstractFloat, Int}, r::Node) =                      Node(2, false, convert(CONST_TYPE, 0f0),       0,      op,  Node(l),       r)
"""
    Node(op::Int, l::Node, r::Union{AbstractFloat, Int})

Short-form to create a scalar/variable node, and apply a binary operator
"""
Node(op::Int, l::Node, r::Union{AbstractFloat, Int}) =                      Node(2, false, convert(CONST_TYPE, 0f0),       0,      op,        l, Node(r))
"""
    Node(op::Int, l::Union{AbstractFloat, Int}, r::Union{AbstractFloat, Int})

Short-form for creating two scalar/variable node, and applying a binary operator
"""
Node(op::Int, l::Union{AbstractFloat, Int}, r::Union{AbstractFloat, Int}) = Node(2, false, convert(CONST_TYPE, 0f0),       0,      op,  Node(l), Node(r))
"""
    Node(val::AbstractFloat)

Create a scalar constant node
"""
Node(val::AbstractFloat) =                                                  Node(convert(CONST_TYPE, val))
"""
    Node(var_string::String)

Create a variable node, using the format `"x1"` to mean feature 1
"""
Node(var_string::String) =                                                  Node(parse(Int, var_string[2:end]))
"""
    Node(var_string::String, varMap::Array{String, 1})

Create a variable node, using a user-passed format
"""
Node(var_string::String, varMap::Array{String, 1}) =                        Node([i for (i, _variable) in enumerate(varMap) if _variable==var_string][1]::Int)


# Copy an equation (faster than deepcopy)
function copyNode(tree::Node)::Node
   if tree.degree == 0
       if tree.constant
           return Node(copy(tree.val))
        else
           return Node(copy(tree.feature))
        end
   elseif tree.degree == 1
       return Node(copy(tree.op), copyNode(tree.l))
    else
        return Node(copy(tree.op), copyNode(tree.l), copyNode(tree.r))
   end
end

function stringOp(op::F, tree::Node, options::Options;
                  bracketed::Bool=false,
                  varMap::Union{Array{String, 1}, Nothing}=nothing)::String where {F}
    if op in [+, -, *, /, ^]
        l = stringTree(tree.l, options, bracketed=false, varMap=varMap)
        r = stringTree(tree.r, options, bracketed=false, varMap=varMap)
        if bracketed
            return "$l $(string(op)) $r"
        else
            return "($l $(string(op)) $r)"
        end
    else
        l = stringTree(tree.l, options, bracketed=true, varMap=varMap)
        r = stringTree(tree.r, options, bracketed=true, varMap=varMap)
        return "$(string(op))($l, $r)"
    end
end

"""
    stringTree(tree::Node, options::Options; kws...)

Convert an equation to a string.

# Arguments

- `varMap::Union{Array{String, 1}, Nothing}=nothing`: what variables
    to print for each feature.
"""
function stringTree(tree::Node, options::Options;
                    bracketed::Bool=false,
                    varMap::Union{Array{String, 1}, Nothing}=nothing)::String
    if tree.degree == 0
        if tree.constant
            return string(tree.val)
        else
            if varMap === nothing
                return "x$(tree.feature)"
            else
                return varMap[tree.feature]
            end
        end
    elseif tree.degree == 1
        return "$(options.unaops[tree.op])($(stringTree(tree.l, options, bracketed=true, varMap=varMap)))"
    else
        return stringOp(options.binops[tree.op], tree, options, bracketed=bracketed, varMap=varMap)
    end
end

# Print an equation
function printTree(io::IO, tree::Node, options::Options; varMap::Union{Array{String, 1}, Nothing}=nothing)
    println(io, stringTree(tree, options, varMap=varMap))
end

function printTree(tree::Node, options::Options; varMap::Union{Array{String, 1}, Nothing}=nothing)
    println(stringTree(tree, options, varMap=varMap))
end

end
