# Partial Equilibrium Model of New and Used Cars
#
# Nested CES demand system with a non-car good:
# - Top level: CES between car aggregate and non-car consumption
# - Car nest: CES between new and used cars
# - New-car nest: CES across brands, then across fuel types within each brand
# - Used-car nest: CES across brands, then across fuel types with habit formation
# - Stock accumulation: used cars = surviving used + depreciated new
#
# Total consumption C[t] and the non-car price p_nc[t] are exogenous.
# The model determines the car/non-car split, allocation across
# new/used, brands, and fuel types, and the market-clearing prices.

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
const b = [:brand1, :brand2, :brand3, :brand4, :brand5]
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
	car_new_b[b,t], "New cars by brand"
	car_new_bf[b,f,t], "New cars by brand and fuel type"
	car_used_b[b,t], "Used-car aggregate by brand"
	car_used_bf[b,f,t], "Used-car stock by brand and fuel type"

	p_C[t], "Price of total consumption"
	p_nc[t], "Price of non-car consumption (exogenous)"
	p_car[t], "User cost of car-service aggregate"
	p_uc_new[t], "User cost of new-car aggregate"
	p_uc_used[t], "User cost of used-car aggregate"
	p_uc_new_b[b,t], "User cost of new car by brand"
	p_uc_used_b[b,t], "User cost of used car by brand"
	p_new_bf[b,f,t], "Purchase price of new car by brand and fuel type (exogenous)"
	p_used_bf[b,f,t], "Price of used car by brand and fuel type"
	p_uc_new_bf[b,f,t], "User cost of new car by brand and fuel type"
	p_uc_used_bf[b,f,t], "User cost of used car by brand and fuel type"

	r[t], "Interest rate"
	δ_new_f[f,t], "First-period depreciation (new→used)"
	δ_used_f[f,t], "Ongoing depreciation of used cars"

	σ_C, "Elasticity: cars vs non-car"
	σ, "Elasticity: new vs used"
	σ_brand, "Elasticity: across brands"
	σ_new, "Elasticity: across fuel types within brand (new)"
	σ_used, "Elasticity: across fuel types within brand (used)"
	μ_car, "CES share: car aggregate"
	μ_noncar, "CES share: non-car consumption"
	μ_new, "CES share: new cars"
	μ_used, "CES share: used cars"
	μ_new_b[b], "CES share: new-car brand"
	μ_used_b[b], "CES share: used-car brand"
	μ_new_f[f], "CES share: fuel type within brand (new)"
	μ_used_f[f], "CES share: fuel type within brand (used)"
	h_f[f], "Habit parameter by fuel type"
	β_h, "Discount on habit premium in user cost (1 = fully forward-looking, 0 = myopic)"

	τ_bf[b,f,t], "Ad-valorem tax (+) or subsidy (-) on new cars by brand and fuel type"
end

# ==============================================================================
# Equations
# ==============================================================================
equations() = @block db begin
	# Top level: cars vs non-car
	car[t = t₁:T],
	car[t] == μ_car * C[t] * (p_car[t] / p_C[t])^(-σ_C)

	C_noncar[t = t₁:T],
	C_noncar[t] == μ_noncar * C[t] * (p_nc[t] / p_C[t])^(-σ_C)

	p_C[t = t₁:T],
	p_C[t] * C[t] == p_car[t] * car[t] + p_nc[t] * C_noncar[t]

	# Car nest: new vs used
	car_new[t = t₁:T],
	car_new[t] == μ_new * car[t] * (p_uc_new[t] / p_car[t])^(-σ)

	car_used[t = t₁:T],
	car_used[t] == μ_used * car[t] * (p_uc_used[t] / p_car[t])^(-σ)

	p_car[t = t₁:T],
	p_car[t] * car[t] == p_uc_new[t] * car_new[t] + p_uc_used[t] * car_used[t]

	# New-car brand nest
	car_new_b[b = b, t = t₁:T],
	car_new_b[b,t] == μ_new_b[b] * car_new[t] * (p_uc_new_b[b,t] / p_uc_new[t])^(-σ_brand)

	p_uc_new[t = t₁:T],
	p_uc_new[t] * car_new[t] == ∑(p_uc_new_b[b,t] * car_new_b[b,t] for b ∈ b)

	# Fuel-type nest within each brand (new)
	car_new_bf[b = b, f = f, t = t₁:T],
	car_new_bf[b,f,t] == μ_new_f[f] * car_new_b[b,t] * (p_uc_new_bf[b,f,t] / p_uc_new_b[b,t])^(-σ_new)

	p_uc_new_b[b = b, t = t₁:T],
	p_uc_new_b[b,t] * car_new_b[b,t] == ∑(p_uc_new_bf[b,f,t] * car_new_bf[b,f,t] for f ∈ f)

	# User cost of new car (b,f): purchase price minus resale value of used car (b,f)
	p_uc_new_bf[b = b, f = f, t = t₁:T-1],
	p_uc_new_bf[b,f,t] == p_new_bf[b,f,t] * (1 + τ_bf[b,f,t]) - (1 - δ_new_f[f,t+1]) / (1 + r[t+1]) * p_used_bf[b,f,t+1]

	p_uc_new_bf[b = b, f = f, t = [T]],
	p_uc_new_bf[b,f,t] == p_new_bf[b,f,t] * (1 + τ_bf[b,f,t]) - (1 - δ_new_f[f,t]) / (1 + r[t]) * p_used_bf[b,f,t]

	# Used-car brand nest
	car_used_b[b = b, t = t₁:T],
	car_used_b[b,t] == μ_used_b[b] * car_used[t] * (p_uc_used_b[b,t] / p_uc_used[t])^(-σ_brand)

	p_uc_used[t = t₁:T],
	p_uc_used[t] * car_used[t] == ∑(p_uc_used_b[b,t] * car_used_b[b,t] for b ∈ b)

	# Fuel-type nest within each brand (used, with habits)
	p_used_bf[b = b, f = f, t = t₁:T],
	car_used_bf[b,f,t] - h_f[f] * car_used_bf[b,f,t-1] == μ_used_f[f] * car_used_b[b,t] * (p_uc_used_bf[b,f,t] / p_uc_used_b[b,t])^(-σ_used)

	p_uc_used_b[b = b, t = t₁:T],
	p_uc_used_b[b,t] * car_used_b[b,t] == ∑(p_uc_used_bf[b,f,t] * (car_used_bf[b,f,t] - h_f[f] * car_used_bf[b,f,t-1]) for f ∈ f)

	p_uc_used_bf[b = b, f = f, t = t₁:T-1],
	p_uc_used_bf[b,f,t] == p_used_bf[b,f,t] - (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * p_used_bf[b,f,t+1] + β_h * (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * h_f[f] * p_uc_used_b[b,t+1] * μ_used_f[f] * (car_used_b[b,t+1] / (car_used_bf[b,f,t+1] - h_f[f] * car_used_bf[b,f,t]))^(1/σ_used)

	p_uc_used_bf[b = b, f = f, t = [T]],
	p_uc_used_bf[b,f,t] == p_used_bf[b,f,t] - (1 - δ_used_f[f,t]) / (1 + r[t]) * p_used_bf[b,f,t] + β_h * (1 - δ_used_f[f,t]) / (1 + r[t]) * h_f[f] * p_uc_used_b[b,t] * μ_used_f[f] * (car_used_b[b,t] / (car_used_bf[b,f,t] - h_f[f] * car_used_bf[b,f,t-1]))^(1/σ_used)

	# Stock accumulation at (b,f) level
	car_used_bf[b = b, f = f, t = t₁:T],
	car_used_bf[b,f,t] == (1 - δ_used_f[f,t]) * car_used_bf[b,f,t-1] + (1 - δ_new_f[f,t]) * car_new_bf[b,f,t-1]
end

steady_state() = @block db begin
	car_new_bf[b = b, f = f, t = [t₀]],
	car_new_bf[b,f,t] == car_new_bf[b,f,t₁]

	car_used_bf[b = b, f = f, t = [t₀]],
	car_used_bf[b,f,t] == car_used_bf[b,f,t₁]
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
optionally extends to the full horizon and re-solves.

Brands are symmetric: each brand gets an equal share of total new and used cars."""
function calibrate(; σ_C_val=0.5, σ_val=3.0, σ_brand_val=5.0, σ_new_val=3.0,
                     σ_used_val=3.0, h_val=0.8, β_h_val=1.0, full_horizon=true)
	data = ModelDictionary(db.model)
	data[σ_C] = σ_C_val
	data[σ] = σ_val
	data[σ_brand] = σ_brand_val
	data[σ_new] = σ_new_val
	data[σ_used] = σ_used_val
	data[h_f] = [h_val, h_val]
	data[β_h] = β_h_val
	data[r] .= 0.04
	data[δ_new_f] .= 0.25
	data[δ_used_f] .= 0.10
	data[C] .= 1.0
	data[p_nc] .= 1.0
	data[p_new_bf] .= 1.0
	data[τ_bf] .= 0.0

	δ₀_cal = 0.25
	δ_cal  = 0.10
	car_share = 0.5

	n_brands = length(b)
	n_fuels = length(f)
	n_bf = n_brands * n_fuels

	# Stock-consistent targets: split car_share equally across (b,f) cells
	car_new_bf_val  = car_share / (n_bf * (1 + (1 - δ₀_cal) / δ_cal * (1 - h_val)))
	car_used_bf_val = (1 - δ₀_cal) / δ_cal * car_new_bf_val
	car_used_bf_adj = car_used_bf_val * (1 - h_val)

	data[car_new_bf[:, :, t₁]] .= car_new_bf_val
	data[car_new_b[:, t₁]] .= car_new_bf_val * n_fuels
	data[car_new[t₁]] = car_new_bf_val * n_bf
	data[car_used_b[:, t₁]] .= car_used_bf_adj * n_fuels
	data[car_used[t₁]] = car_used_bf_adj * n_bf
	data[car[t₁]] = car_share
	data[C_noncar[t₁]] = data[C[t₁]] - car_share
	data[p_used_bf[:, :, t₁]] .= 0.5

	global T = t₁
	cal = equations() + steady_state()
	# Swap μ_new_f/μ_used_f against a single brand (symmetric → same result for all)
	@endo_exo_swap! cal begin
		μ_car, car[t₁]
		μ_noncar, C_noncar[t₁]
		μ_new_b, car_new_b[b, t₁]
		μ_used_b, car_used_b[b, t₁]
		μ_new_f, car_new_bf[:brand1, f, t₁]
		μ_used_f, p_used_bf[:brand1, f, t₁]
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
	println("  μ_new_b = ", [result[μ_new_b[bb]] for bb in b])
	println("  μ_used_b = ", [result[μ_used_b[bb]] for bb in b])
	println("  μ_new_f = ", [result[μ_new_f[ff]] for ff in f])
	println("  μ_used_f = ", [result[μ_used_f[ff]] for ff in f])

	return result
end
