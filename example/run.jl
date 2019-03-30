
using MathOptFormat
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

DIR_PATH = dirname(@__FILE__)

push!(Base.LOAD_PATH, joinpath(DIR_PATH, "..",".."))

using ModelCompare

cd(DIR_PATH)

model = MathOptFormat.MOF.Model()
x = MOI.add_variable(model)
MOI.set(model, MOI.VariableName(), x, "x")
c = MOI.add_constraint(model,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}.([1],[x]), 0.0),
    MOI.GreaterThan(1.0)
)
MOI.set(model, MOI.ConstraintName(), c, "c")
c2 = MOI.add_constraint(model,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
    MOI.Interval(-1.0, 1.0)
)
MOI.set(model, MOI.ConstraintName(), c2, "c2")
MOI.write_to_file(model, "orig_A.mps")

remove_intervals(model)
MOI.write_to_file(model, "orig_A_no_interval.mps")

# ------------------------------------------------------------------

model_ = MathOptFormat.MOF.Model()
z = MOI.add_variable(model_)
MOI.set(model_, MOI.VariableName(), z, "z")
x = MOI.add_variable(model_)
MOI.set(model_, MOI.VariableName(), x, "x")
c = MOI.add_constraint(model_,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}.([1,2],[z,x]), 0.0),
    MOI.GreaterThan(1.0)
)
MOI.set(model_, MOI.ConstraintName(), c, "c")
c2 = MOI.add_constraint(model_,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
    MOI.Interval(-1.0, 1.0)
)
MOI.set(model_, MOI.ConstraintName(), c2, "c2")
MOI.write_to_file(model_, "orig_B.mps")
remove_intervals(model_)



compare_models(model, model_)