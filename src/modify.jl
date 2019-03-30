ctr_name(::MOI.EqualTo) = "FX"
ctr_name(::MOI.LessThan) = "UB"
ctr_name(::MOI.GreaterThan) = "LB"
ctr_name(::MOI.Integer) = "IN"
ctr_name(::MOI.ZeroOne) = "BN"

function fill_constraint_names(model)
    for (F,S) in MOI.get(model, MOI.ListOfConstraints())
        if F <: MOI.ScalarAffineFunction#, MOI.ScalarQuadraticFunction]
            # check
        elseif F <: MOI.VectorOfVariables
            # error
        elseif F <: MOI.SingleVariable
            for index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
                name = MOI.get(model, MOI.ConstraintName(), index)
                if name == ""
                    f = MOI.get(model, MOI.ConstraintFunction(), index)
                    name = MOI.get(model, MOI.VariableName(), f.variable)
                    MOI.set(model, MOI.ConstraintName(), index, name*"_"*ctr_name(S))
                end
            end
        else
            #error
        end
    end
end

function remove_intervals(model)
    for (F,S) in MOI.get(model, MOI.ListOfConstraints())
        if S <: MOI.Interval
            for index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
                cname = MOI.get(model, MOI.ConstraintName(), index)
                func = MOI.get(model, MOI.ConstraintFunction(), index)
                set = MOI.get(model, MOI.ConstraintSet(), index)
                c0 = MOI.transform(model, index, MOI.LessThan(set.upper))
                MOI.set(model, MOI.ConstraintName(), c0, cname) # choose a nice name
                c = MOI.add_constraint(model, func, MOI.GreaterThan(set.lower))
                MOI.set(model, MOI.ConstraintName(), c, cname * "_LB") # choose a nice name
            end
        end
    end
end