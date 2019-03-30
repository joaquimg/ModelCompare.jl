function sort_by_name(model::MOI.ModelLike, func::MOI.AbstractFunction)
    error("`sort_by_name` is not defined for $(typeof(func))")
end
function sort_by_name(model::MOI.ModelLike, func::MOI.SingleVariable)
end
function sort_by_name(model::MOI.ModelLike, func::MOI.ScalarAffineFunction)
    names = Array{String}(undef, length(func.terms))
    for i in 1:length(func.terms)
        names[i] = MOI.get(model, MOI.VariableName(), func.terms[i].variable_index)
    end
    perm = sortperm(names)
    func.terms = func.terms[perm]
    return
end
function sort_by_name(model::MOI.ModelLike, ce::Vector{CE}) where CE <:MOIU.ConstraintEntry
    names = Array{String}(undef, length(ce))
    for i in 1:length(ce)
        names[i] = MOI.get(model, MOI.ConstraintName(), ce[i][1])
    end
    perm = sortperm(names)
    ce .= ce[perm]
    return
end
function sort_by_variable_name(model::MOI.ModelLike, ce::Vector{CE}) where CE <: MOIU.ConstraintEntry
    names = Array{String}(undef, length(ce))
    for i in 1:length(ce)
        f = MOI.get(model, MOI.ConstraintFunction(), ce[i][1])
        names[i] = MOI.get(model, MOI.VariableName(), f.variable)
    end
    perm = sortperm(names)
    ce .= ce[perm]
    return
end
function sort_model(model_u)
    model = model_u.model
    @show fieldnames(typeof(model))
    sort_by_name(model, model.objective)
    sv = model.moi_singlevariable
    for f in fieldnames(typeof(sv))
        sort_by_variable_name(model, getfield(sv, f))
    end
    saf = model.moi_scalaraffinefunction
    for f in fieldnames(typeof(saf))
        vec = getfield(saf, f)
        sort_by_name(model, vec)
        for ce in vec
            sort_by_name(model,  ce[2])
        end
    end
    return
end