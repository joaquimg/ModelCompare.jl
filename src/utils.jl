
function pass_variables(m1, m2, f::MOI.AbstractFunction, map, tol::Number)
    vis = non_zero_variables(f, tol)
    pass_variables(m1, m2, vis, map)
    return
end

function pass_variables(m1, m2, vis, map)
    for vi in vis
        if !haskey(map, vi)
            vname = MOI.get(m1, MOI.VariableName(), vi)
            x = MOI.add_variable(m2)
            MOI.set(m2, MOI.VariableName(), x, vname)
            map[vi] = x
        end
    end
    return
end

is_approx_zero(lhs::MOI.EqualTo, tol = 1e-8) = abs(lhs.value) < tol
is_approx_zero(lhs::MOI.GreaterThan, tol = 1e-8) = abs(lhs.lower) < tol
is_approx_zero(lhs::MOI.LessThan, tol = 1e-8) = abs(lhs.upper) < tol
is_approx_zero(lhs::S, tol = 1e-8) where S<:MOI.AbstractSet = true

set_difference(lhs::S, rhs::S) where S = lhs
set_sum(lhs::S, rhs::S) where S = lhs

set_sum(lhs::MOI.EqualTo, rhs::MOI.EqualTo) = MOI.EqualTo(rhs.value + lhs.value)
set_difference(lhs::MOI.EqualTo, rhs::MOI.EqualTo) = MOI.EqualTo(lhs.value - rhs.value)

set_sum(lhs::MOI.GreaterThan, rhs::MOI.GreaterThan) = MOI.GreaterThan(rhs.lower + lhs.lower)
set_difference(lhs::MOI.GreaterThan, rhs::MOI.GreaterThan) = MOI.GreaterThan(lhs.lower - rhs.lower)

set_sum(lhs::MOI.LessThan, rhs::MOI.LessThan) = MOI.LessThan(rhs.upper + lhs.upper)
set_difference(lhs::MOI.LessThan, rhs::MOI.LessThan) = MOI.LessThan(lhs.upper - rhs.upper)

function non_zero_variables(f::MOI.ScalarAffineFunction, tol)
    vis = Set{VI}()
    to_stay = Int[]
    for i in 1:length(f.terms)
        if abs(MOI.coefficient(f.terms[i])) > tol
            push!(vis, f.terms[i].variable_index)
            push!(to_stay, i)
        end
    end
    f.terms = f.terms[to_stay]
    return vis
end
function non_zero_variables(f::MOI.SingleVariable, tol)
    vis = Set{VI}(f.variable)
    return vis
end

# TODO - Other functions
function unsafe_move_variables(set::Set, f1, f2)
    for i in set
        unsafe_move_variable(i, f1, f2)
    end
end
function unsafe_move_variable(vi, f1::MOI.ScalarAffineFunction, f2::MOI.ScalarAffineFunction) # suppose canonical
    for i in 1:length(f1.terms)
        if MOIU._hasvar(f1.terms[i], vi)
            push!(f2.terms, f1.terms[i])
            deleteat!(f1.terms, i) # suppose canonical
            return true
        end
    end
    return false
end
