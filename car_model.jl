# Partial Equilibrium Model of New and Used Cars
#
# Nested CES demand system with a non-car good:
# - Top level: CES between car aggregate and non-car consumption
# - Car nest: CES between new and used cars
# - New-car nest: CES across fuel types
# - Used-car nest: CES across fuel types with habit formation
# - Stock accumulation: used cars = surviving used + depreciated new
#
# Total consumption C[t] and the non-car price p_nc[t] are exogenous.
# The model determines the car/non-car split, allocation across
# new/used and fuel types, and the market-clearing prices.

import JuMP
using JuMP: Model, set_silent, name
using Ipopt
using SquareModels
using CairoMakie

# ==============================================================================
# Model container and time/set configuration
# ==============================================================================
db = ModelDictionary(Model(Ipopt.Optimizer))
set_silent(db.model)

const f = [:petrol, :electric]
const t₀ = 2024
const max_T = 2100
const t = t₀:max_T

t₁::Int = t₀ + 1
T::Int = max_T

# ==============================================================================
# Variables
# ==============================================================================
@variables db.model begin
	C[t], "Total consumption (exogenous)"
	C_noncar[t], "Non-car consumption"
	car[t], "Car-service aggregate"
	car_new[t], "New-car aggregate"
	car_used[t], "Used-car aggregate"
	car_new_f[f,t], "New cars by fuel type"
	car_used_f[f,t], "Used-car stock by fuel type"

	p_C[t], "Price of total consumption"
	p_nc[t], "Price of non-car consumption (exogenous)"
	p_car[t], "User cost of car-service aggregate"
	p_uc_new[t], "User cost of new-car aggregate"
	p_uc_used[t], "User cost of used-car aggregate"
	p_new_f[f,t], "Purchase price of new car by fuel type (exogenous)"
	p_used_f[f,t], "Price of used car by fuel type"
	p_uc_new_f[f,t], "User cost of new car by fuel type"
	p_uc_used_f[f,t], "User cost of used car by fuel type"

	r[t], "Interest rate"
	δ_new_f[f,t], "First-period depreciation (new→used)"
	δ_used_f[f,t], "Ongoing depreciation of used cars"

	σ_C, "Elasticity: cars vs non-car"
	σ, "Elasticity: new vs used"
	σ_new, "Elasticity: across fuel types (new)"
	σ_used, "Elasticity: across fuel types (used)"
	μ_car, "CES share: car aggregate"
	μ_noncar, "CES share: non-car consumption"
	μ_new, "CES share: new cars"
	μ_used, "CES share: used cars"
	μ_new_f[f], "CES share: new car by fuel type"
	μ_used_f[f], "CES share: used car by fuel type"
	h_f[f], "Habit parameter by fuel type"
	β_h, "Discount on habit premium in user cost (1 = fully forward-looking, 0 = myopic)"

	τ_f[f,t], "Ad-valorem tax (+) or subsidy (-) on new cars by fuel type"
end

# ==============================================================================
# Equations
# ==============================================================================
equations() = @block db begin
	car[t = t₁:T],
	car[t] == μ_car * C[t] * (p_car[t] / p_C[t])^(-σ_C)

	C_noncar[t = t₁:T],
	C_noncar[t] == μ_noncar * C[t] * (p_nc[t] / p_C[t])^(-σ_C)

	p_C[t = t₁:T],
	p_C[t] * C[t] == p_car[t] * car[t] + p_nc[t] * C_noncar[t]

	car_new[t = t₁:T],
	car_new[t] == μ_new * car[t] * (p_uc_new[t] / p_car[t])^(-σ)

	car_used[t = t₁:T],
	car_used[t] == μ_used * car[t] * (p_uc_used[t] / p_car[t])^(-σ)

	p_car[t = t₁:T],
	p_car[t] * car[t] == p_uc_new[t] * car_new[t] + p_uc_used[t] * car_used[t]

	car_new_f[f = f, t = t₁:T],
	car_new_f[f,t] == μ_new_f[f] * car_new[t] * (p_uc_new_f[f,t] / p_uc_new[t])^(-σ_new)

	p_uc_new[t = t₁:T],
	p_uc_new[t] * car_new[t] == ∑(p_uc_new_f[f,t] * car_new_f[f,t] for f ∈ f)

	p_uc_new_f[f = f, t = t₁:T-1],
	p_uc_new_f[f,t] == p_new_f[f,t] * (1 + τ_f[f,t]) - (1 - δ_new_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1]

	p_uc_new_f[f = f, t = [T]],
	p_uc_new_f[f,t] == p_new_f[f,t] * (1 + τ_f[f,t]) - (1 - δ_new_f[f,t]) / (1 + r[t]) * p_used_f[f,t]

	car_used_f[f = f, t = t₁:T],
	car_used_f[f,t] == (1 - δ_used_f[f,t]) * car_used_f[f,t-1] + (1 - δ_new_f[f,t]) * car_new_f[f,t-1]

	p_used_f[f = f, t = t₁:T],
	car_used_f[f,t] - h_f[f] * car_used_f[f,t-1] == μ_used_f[f] * car_used[t] * (p_uc_used_f[f,t] / p_uc_used[t])^(-σ_used)

	p_uc_used[t = t₁:T],
	p_uc_used[t] * car_used[t] == ∑(p_uc_used_f[f,t] * (car_used_f[f,t] - h_f[f] * car_used_f[f,t-1]) for f ∈ f)

	p_uc_used_f[f = f, t = t₁:T-1],
	p_uc_used_f[f,t] == p_used_f[f,t] - (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1] + β_h * (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * h_f[f] * p_uc_used[t+1] * μ_used_f[f] * (car_used[t+1] / (car_used_f[f,t+1] - h_f[f] * car_used_f[f,t]))^(1/σ_used)

	p_uc_used_f[f = f, t = [T]],
	p_uc_used_f[f,t] == p_used_f[f,t] - (1 - δ_used_f[f,t]) / (1 + r[t]) * p_used_f[f,t] + β_h * (1 - δ_used_f[f,t]) / (1 + r[t]) * h_f[f] * p_uc_used[t] * μ_used_f[f] * (car_used[t] / (car_used_f[f,t] - h_f[f] * car_used_f[f,t-1]))^(1/σ_used)
end

steady_state() = @block db begin
	car_new_f[f = f, t = [t₀]],
	car_new_f[f,t] == car_new_f[f,t₁]

	car_used_f[f = f, t = [t₀]],
	car_used_f[f,t] == car_used_f[f,t₁]
end

"""Copy endogenous variable values from period `t_source` to each period in `periods`,
using the block's endogenous list to discover all time-indexed variables automatically."""
function extend_to_horizon!(data::ModelDictionary, block::Block, t_source, periods)
	suffix = "$t_source]"
	for var in block.variables
		k = name(var)
		endswith(k, suffix) || continue
		val = data[k]
		val === nothing && continue
		for tt in periods
			data[k[1:end-length(suffix)] * "$tt]"] = val
		end
	end
end

# ==============================================================================
# Calibration
# ==============================================================================
"""Calibrate the model at a steady state and return the solved baseline data.

Keyword arguments override the default parameter values. The calibration solves
for share parameters (μ's) given stock-consistent quantity targets, then
optionally extends to the full horizon and re-solves."""
function calibrate(; σ_C_val=0.5, σ_val=3.0, σ_new_val=3.0, σ_used_val=3.0,
                     h_val=0.8, β_h_val=1.0, full_horizon=true)
	data = ModelDictionary(db.model)
	data[σ_C] = σ_C_val
	data[σ] = σ_val
	data[σ_new] = σ_new_val
	data[σ_used] = σ_used_val
	data[h_f] = [h_val, h_val]
	data[β_h] = β_h_val
	data[r] .= 0.04
	data[δ_new_f] .= 0.25
	data[δ_used_f] .= 0.10
	data[C] .= 1.0
	data[p_nc] .= 1.0
	data[p_new_f] .= 1.0
	data[τ_f] .= 0.0

	δ₀_cal = 0.25
	δ_cal  = 0.10
	car_share = 0.5

	n_fuels = length(f)
	car_new_f_val  = car_share / (n_fuels * (1 + (1 - δ₀_cal) / δ_cal * (1 - h_val)))
	car_used_f_val = (1 - δ₀_cal) / δ_cal * car_new_f_val
	car_used_adj   = car_used_f_val * (1 - h_val)

	data[car_new_f[:, t₁]] .= car_new_f_val
	data[car_new[t₁]] = car_new_f_val * n_fuels
	data[car_used[t₁]] = car_used_adj * n_fuels
	data[car[t₁]] = car_share
	data[C_noncar[t₁]] = data[C[t₁]] - car_share
	data[p_used_f[:, t₁]] .= 0.5

	global T = t₁
	cal = equations() + steady_state()
	@endo_exo_swap! cal begin
		μ_car, car[t₁]
		μ_noncar, C_noncar[t₁]
		μ_new_f, car_new_f[f, t₁]
		μ_used_f, p_used_f[f, t₁]
		μ_new, car_new[t₁]
		μ_used, car_used[t₁]
	end
	result = solve(cal, data; replace_nothing=1.0)

	if full_horizon
		global T = max_T
		extend_to_horizon!(result, cal, t₁, t₁+1:T)
		solve!(equations(), result)
	else
		global T = max_T
	end

	println("Calibrated parameters:")
	println("  μ_car = ", result[μ_car])
	println("  μ_noncar = ", result[μ_noncar])
	println("  μ_new = ", result[μ_new])
	println("  μ_used = ", result[μ_used])
	println("  μ_new_f = ", [result[μ_new_f[ff]] for ff in f])
	println("  μ_used_f = ", [result[μ_used_f[ff]] for ff in f])

	return result
end
