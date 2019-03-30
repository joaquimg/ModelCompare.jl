module ModelCompare

using MathOptFormat
using MathOptInterface

export compare_models, fill_constraint_names, remove_intervals

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

include("sorting.jl")
include("modify.jl")
include("utils.jl")

function compare_models(model1::MOI.ModelLike, model2::MOI.ModelLike, tol = 1e-8)
    io = STDIO

    # TODO: give test-friendly feedback instead of errors?
    extra_1_model = MathOptFormat.MOF.Model()
    extra_2_model = MathOptFormat.MOF.Model()

    model1_vnames = [MOI.get(model, MOI.VariableName(), index) 
                     for index in MOI.get(model, MOI.ListOfVariableIndices())]
    model1_cnames = String[]
    for (F,S) in MOI.get(model, MOI.ListOfConstraints())
        for index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            push!(model1_cnames, MOI.get(model, MOI.ConstraintName(), index))
        end
    end
    # model1_cnames = vec([MOI.get(model, MOI.ConstraintName(), index) 
    #                      for (F,S) in MOI.get(model, MOI.ListOfConstraints()),
    #                          index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())])

    extra_1_vnames, extra_2_vnames, same_vnames = compare_variablenames(model2, model1_vnames)
    println("Variables in both: $(length(same_vnames))")
    println("Variables only in left: $(length(extra_1_vnames))")
    println("Variables only in right: $(length(extra_2_vnames))")
    extra_1_cnames, extra_2_cnames, same_cnames = compare_constraintnames(model2, model1_cnames)
    println("Constraints in both: $(length(same_cnames))")
    println("Constraints only in left: $(length(extra_1_cnames))")
    println("Constraints only in right: $(length(extra_2_cnames))")

    variablemap_2to1 = Dict{VI,VI}()
    variablemap_1to2 = Dict{VI,VI}()
    for vname in same_vnames
        index1 = MOI.get(model1, VI, vname)
        index2 = MOI.get(model2, VI, vname)
        variablemap_2to1[index2] = index1
        variablemap_1to2[index1] = index2
    end

    variablemap_1_old_to_new = Dict{VI,VI}()
    extra_1_vi = Set{VI}()
    for vname in extra_1_vnames
        index1 = MOI.get(model1, VI, vname)
        # add variables
        x = MOI.add_variable(extra_1_model)
        MOI.set(extra_1_model, MOI.VariableName(), x, vname)
        variablemap_1_old_to_new[index1] = x
        push!(extra_1_vi, index1)
    end
    
    variablemap_2_old_to_new = Dict{VI,VI}()
    extra_2_vi = Set{VI}()
    for vname in extra_2_vnames
        index2 = MOI.get(model2, VI, vname)
        # add variables
        x = MOI.add_variable(extra_2_model)
        MOI.set(extra_2_model, MOI.VariableName(), x, vname)
        variablemap_2_old_to_new[index2] = x
        push!(extra_2_vi, index2)
    end

    for cname in same_cnames
        index1 = MOI.get(model1, CI, cname)
        index2 = MOI.get(model2, CI, cname)
        f1 = MOI.get(model1, MOI.ConstraintFunction(), index1)
        f2 = MOI.get(model2, MOI.ConstraintFunction(), index2)

        f1_type = typeof(f1)
        f2_type = typeof(f2)
        @assert f1_type == f2_type

        s1 = MOI.get(model1, MOI.ConstraintSet(), index1)
        s2 = MOI.get(model2, MOI.ConstraintSet(), index2)
        s_12 = set_difference(s1,s2)
        s_21 = set_difference(s2,s1)
        is_approx_zero_s_12 = is_approx_zero(s_12, tol)

        if f1_type <: MOI.ScalarAffineFunction#, MOI.ScalarQuadraticFunction]

            f1_extra = MOI.zero(f1_type)
            f2_extra = MOI.zero(f2_type)
            unsafe_move_variables(extra_1_vi, f1, f1_extra)
            unsafe_move_variables(extra_2_vi, f2, f2_extra)

            f_12 = MOIU.canonical(f1 - MOIU.mapvariables(variablemap_2to1, f2))
            f_21 = MOI.zero(f1_type)

            if !MOIU.isapprox_zero(f_12, tol)
                f_21 = MOIU.canonical(f2 - MOIU.mapvariables(variablemap_1to2, f1))
                pass_variables(model1, extra_1_model, f_12, variablemap_1_old_to_new, tol)
                pass_variables(model2, extra_2_model, f_21, variablemap_2_old_to_new, tol)
                f_12 = MOIU.isapprox_zero(f1_extra, tol) ? f_12 : MOIU.canonical(f_12 + f1_extra)
                f_21 = MOIU.isapprox_zero(f2_extra, tol) ? f_21 : MOIU.canonical(f_21 + f2_extra)
                c1 = MOI.add_constraint(extra_1_model, MOIU.mapvariables(variablemap_1_old_to_new, f_12), s_12)
                c2 = MOI.add_constraint(extra_2_model, MOIU.mapvariables(variablemap_2_old_to_new, f_21), s_21)
                MOI.set(extra_1_model, MOI.ConstraintName(), c1, cname)
                MOI.set(extra_2_model, MOI.ConstraintName(), c2, cname)
            else
                if !MOIU.isapprox_zero(f1_extra, tol)
                    c1 = MOI.add_constraint(extra_1_model, MOIU.mapvariables(variablemap_1_old_to_new, f1_extra), s_12)
                    MOI.set(extra_1_model, MOI.ConstraintName(), c1, cname)
                end
                if !MOIU.isapprox_zero(f2_extra, tol)
                    c2 = MOI.add_constraint(extra_2_model, MOIU.mapvariables(variablemap_2_old_to_new, f2_extra), s_21)
                    MOI.set(extra_2_model, MOI.ConstraintName(), c2, cname)
                end
            end
        elseif f1_type <: MOI.SingleVariable
            if f1.variable != MOIU.mapvariables(variablemap_2to1, f2).variable || !is_approx_zero_s_12
                pass_variables(model1, extra_1_model, f1, variablemap_1_old_to_new, tol)
                pass_variables(model2, extra_2_model, f2, variablemap_2_old_to_new, tol)
                c1 = MOI.add_constraint(extra_1_model, f1_extra, s_12)
                c2 = MOI.add_constraint(extra_2_model, f2_extra, s_21)
                MOI.set(extra_1_model, MOI.ConstraintName(), c1, cname)
                MOI.set(extra_2_model, MOI.ConstraintName(), c2, cname)
            end
        elseif f1_type == MOI.VectorOfVariables
            error("Vector os variables comparisson not supported")
        else
            error("Invalid function type")
        end
    end
    for cname in extra_1_cnames
        index = MOI.get(model1, CI, cname)
        s = MOI.get(model1, MOI.ConstraintSet(), index)
        f = MOI.get(model1, MOI.ConstraintFunction(), index)
        pass_variables(model1, extra_1_model, f, variablemap_1_old_to_new, tol)
        c = MOI.add_constraint(extra_1_model, MOIU.mapvariables(variablemap_1_old_to_new, f), s)
        MOI.set(extra_2_model, MOI.ConstraintName(), c, cname)
    end
    for cname in extra_2_cnames
        index = MOI.get(model2, CI, cname)
        s = MOI.get(model2, MOI.ConstraintSet(), index)
        f = MOI.get(model2, MOI.ConstraintFunction(), index)
        pass_variables(model2, extra_2_model, f, variablemap_2_old_to_new, tol)
        c = MOI.add_constraint(extra_2_model, MOIU.mapvariables(variablemap_2_old_to_new, f), s)
        MOI.set(extra_2_model, MOI.ConstraintName(), c, cname)
    end

    sort_model(extra_1_model)
    sort_model(extra_2_model)

    if true
        lp_1 = MathOptFormat.LP.Model()
        MOI.copy_to(lp_1, extra_1_model)
        lp_2 = MathOptFormat.LP.Model()
        MOI.copy_to(lp_2, extra_2_model)
        MOI.write_to_file(lp_1, "A.lp")
        MOI.write_to_file(lp_2, "B.lp")
    end

    if true
        mps_1 = MathOptFormat.MPS.Model()
        MOI.copy_to(mps_1, extra_1_model)
        mps_2 = MathOptFormat.MPS.Model()
        MOI.copy_to(mps_2, extra_2_model)
        MOI.write_to_file(mps_1, "A.mps")
        MOI.write_to_file(mps_2, "B.mps")
    end

    MOI.write_to_file(extra_1_model, "A.mof")
    MOI.write_to_file(extra_2_model, "B.mof")

end


function compare_variablenames(model, variablenames)
    missing_names = String[] # TODO use sets?
    extra_names = String[]
    same_names = String[]
    seen_name = Dict(name => false for name in variablenames)
    for index in MOI.get(model, MOI.ListOfVariableIndices())
        vname = MOI.get(model, MOI.VariableName(), index)
        if !haskey(seen_name, vname)
            push!(extra_names, vname)
            # error("Variable with name $vname present in model but not expected list of variable names.")
        elseif seen_name[vname]
            println("Variable with name $vname present twice in model (shouldn't happen!)")
        else
            seen_name[vname] = true
            push!(same_names, vname)
        end
    end
    for (vname,seen) in seen_name
        if !seen
            push!(missing_names, vname)
            # error("Did not find variable with name $vname in intance.")
        end
    end
    return missing_names, extra_names, same_names
end
function compare_constraintnames(model, constraintnames)
    missing_names = String[]
    extra_names = String[]
    same_names = String[]
    seen_name = Dict(name => false for name in constraintnames)
    for (F,S) in MOI.get(model, MOI.ListOfConstraints())
        for index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            cname = MOI.get(model, MOI.ConstraintName(), index)
            if !haskey(seen_name, cname)
                push!(extra_names, cname)
                # error("Constraint with name $cname present in model but not expected list of constraint names.")
            elseif seen_name[cname]
                println("Constraint with name $cname present twice in model (shouldn't happen!)")
            else
                seen_name[cname] = true
                push!(same_names, cname)
            end
        end
    end
    for (cname,seen) in seen_name
        if !seen
            push!(missing_names, cname)
            # error("Did not find constraint with name $cname in intance.")
        end
    end
    return missing_names, extra_names, same_names
end


end
