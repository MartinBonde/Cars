# Partial Equilibrium Model of New and Used Cars
#
# Implements the CES demand system from cars.lyx with a non-car good:
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
using JuMP: Model, set_silent
using Ipopt
using SquareModels
using CairoMakie

# ==============================================================================
# Model container and time/set configuration
# ==============================================================================
db = ModelDictionary(Model(Ipopt.Optimizer))
set_silent(db.model)

const f = [:petrol, :electric]  # Fuel types
const t₀ = 2024
const max_T = 2060
const t = t₀:max_T

t₁::Int = t₀ + 1
T::Int = max_T

# ==============================================================================
# Variables
# ==============================================================================
@variables db.model begin
	# Quantities — top level
	C[t], "Total consumption (exogenous)"
	c_nc[t], "Non-car consumption"
	d[t], "Car-service aggregate"
	d_new[t], "New-car aggregate"
	d_used[t], "Used-car aggregate"
	d_new_f[f,t], "New cars by fuel type"
	d_used_f[f,t], "Used-car stock by fuel type"

	# Prices
	p_C[t], "Price of total consumption"
	p_nc[t], "Price of non-car consumption (exogenous)"
	p_d[t], "Shadow price of car-service aggregate"
	p_uc_new[t], "User cost of new-car aggregate"
	p_uc_used[t], "User cost of used-car aggregate"
	p_new_f[f,t], "Purchase price of new car by fuel type (exogenous)"
	p_used_f[f,t], "Price of used car by fuel type"
	p_uc_new_f[f,t], "User cost of new car by fuel type"
	p_uc_f[f,t], "User cost of used car by fuel type"
	p_uc_d[t], "User cost of car-service aggregate"

	# Exogenous
	r[t], "Interest rate"
	δ_d[t], "Depreciation of car-service aggregate"
	δ0_f[f,t], "First-period depreciation (new→used)"
	δ_f[f,t], "Ongoing depreciation of used cars"

	# Parameters / calibrated
	σ_C, "Elasticity: cars vs non-car"
	σ, "Elasticity: new vs used"
	σ_new, "Elasticity: across fuel types (new)"
	σ_used, "Elasticity: across fuel types (used)"
	μ_d, "CES share: car aggregate"
	μ_nc, "CES share: non-car consumption"
	μ_new, "CES share: new cars"
	μ_used, "CES share: used cars"
	μ_new_f[f], "CES share: new car by fuel type"
	μ_used_f[f], "CES share: used car by fuel type"
	h_f[f], "Habit parameter by fuel type"
end

# ==============================================================================
# Equations
# ==============================================================================
equations = @block db begin
	# ------------------------------------------------------------------
	# Top level: car aggregate vs non-car consumption
	# ------------------------------------------------------------------
	d[t = t₁:T],
	d[t] == μ_d * C[t] * (p_uc_d[t] / p_C[t])^(-σ_C)

	c_nc[t = t₁:T],
	c_nc[t] == μ_nc * C[t] * (p_nc[t] / p_C[t])^(-σ_C)

	p_C[t = t₁:T],
	p_C[t] * C[t] == p_uc_d[t] * d[t] + p_nc[t] * c_nc[t]

	# ------------------------------------------------------------------
	# Car nest: new vs used
	# ------------------------------------------------------------------
	d_new[t = t₁:T],
	d_new[t] == μ_new * d[t] * (p_uc_new[t] / p_d[t])^(-σ)

	d_used[t = t₁:T],
	d_used[t] == μ_used * d[t] * (p_uc_used[t] / p_d[t])^(-σ)

	p_d[t = t₁:T],
	p_d[t] * d[t] == p_uc_new[t] * d_new[t] + p_uc_used[t] * d_used[t]

	p_uc_d[t = t₁:T-1],
	p_uc_d[t] == p_d[t] - (1 - δ_d[t+1]) / (1 + r[t+1]) * p_d[t+1]

	p_uc_d[t = [T]],
	p_uc_d[t] == p_d[t] * (1 - (1 - δ_d[t]) / (1 + r[t]))

	# ------------------------------------------------------------------
	# New-car nest: across fuel types
	# ------------------------------------------------------------------
	d_new_f[f = f, t = t₁:T],
	d_new_f[f,t] == μ_new_f[f] * d_new[t] * (p_uc_new_f[f,t] / p_uc_new[t])^(-σ_new)

	p_uc_new[t = t₁:T],
	p_uc_new[t] * d_new[t] == ∑(p_uc_new_f[f,t] * d_new_f[f,t] for f ∈ f)

	p_uc_new_f[f = f, t = t₁:T-1],
	p_uc_new_f[f,t] == p_new_f[f,t] - (1 - δ0_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1]

	p_uc_new_f[f = f, t = [T]],
	p_uc_new_f[f,t] == p_new_f[f,t] - (1 - δ0_f[f,t]) / (1 + r[t]) * p_used_f[f,t]

	# ------------------------------------------------------------------
	# Used-car nest: across fuel types (with habits)
	# ------------------------------------------------------------------
	d_used_f[f = f, t = t₁:T],
	d_used_f[f,t] == (1 - δ_f[f,t]) * d_used_f[f,t-1] + (1 - δ0_f[f,t]) * d_new_f[f,t-1]

	# Price adjusts to clear the used-car market given stocks
	p_used_f[f = f, t = t₁:T],
	d_used_f[f,t] - h_f[f] * d_used_f[f,t-1] == μ_used_f[f] * d_used[t] * (p_uc_f[f,t] / p_uc_used[t])^(-σ_used)

	p_uc_used[t = t₁:T],
	p_uc_used[t] * d_used[t] == ∑(p_uc_f[f,t] * (d_used_f[f,t] - h_f[f] * d_used_f[f,t-1]) for f ∈ f)

	p_uc_f[f = f, t = t₁:T-1],
	p_uc_f[f,t] == p_used_f[f,t] - (1 - δ_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1] + (1 - δ_f[f,t+1]) / (1 + r[t+1]) * h_f[f] * p_uc_used[t+1] * μ_used_f[f] * (d_used[t+1] / (d_used_f[f,t+1] - h_f[f] * d_used_f[f,t]))^(1/σ_used)

	p_uc_f[f = f, t = [T]],
	p_uc_f[f,t] == p_used_f[f,t] - (1 - δ_f[f,t]) / (1 + r[t]) * p_used_f[f,t] + (1 - δ_f[f,t]) / (1 + r[t]) * h_f[f] * p_uc_used[t] * μ_used_f[f] * (d_used[t] / (d_used_f[f,t] - h_f[f] * d_used_f[f,t-1]))^(1/σ_used)
end

# ==============================================================================
# Data
# ==============================================================================
# Elasticities
db[σ_C] = 0.5
db[σ] = 1.5
db[σ_new] = 3.0
db[σ_used] = 3.0

# Share parameters (will be calibrated)
db[μ_d] = 0.2
db[μ_nc] = 0.8
db[μ_new] = 0.3
db[μ_used] = 0.7
db[μ_new_f] = [0.5, 0.5]   # petrol, electric
db[μ_used_f] = [0.5, 0.5]

# Habits: h = 1 - 1/avg_holding_period, ~8 years average
db[h_f] = [0.875, 0.875]

# Interest rate and depreciation
db[r] .= 0.03
db[δ_d] .= 0.10
db[δ0_f] .= 0.25    # New cars lose 25% in first year
db[δ_f] .= 0.10     # Used cars depreciate 10% per year

# Exogenous
db[C] .= 1000.0
db[p_nc] .= 1.0
db[p_new_f] .= 1.0

# ==============================================================================
# Starting values (hot-start for solver)
# ==============================================================================
# In steady state: d_used_f = (1-δ0) / δ * d_new_f = 0.75/0.10 * d_new_f = 7.5 * d_new_f
db[d_new_f[:petrol, t]] .= 10.0
db[d_new_f[:electric, t]] .= 10.0
db[d_used_f[:petrol, t]] .= 75.0
db[d_used_f[:electric, t]] .= 75.0

db[d_new] .= 20.0
db[d_used] .= 150.0
db[d] .= 170.0
db[c_nc] .= 830.0

db[p_used_f] .= 0.5
db[p_d] .= 0.1
db[p_C] .= 1.0
db[p_uc_new_f] .= 0.64   # ≈ 1 - 0.75/1.03 * 0.5
db[p_uc_new] .= 0.64
db[p_uc_f] .= 0.06       # ≈ 0.5 - 0.9/1.03 * 0.5
db[p_uc_used] .= 0.06
db[p_uc_d] .= 0.01

# ==============================================================================
# Calibration
# ==============================================================================
calibration = copy(equations)
@endo_exo_swap! calibration begin
	μ_d, d[t₁]
	μ_nc, c_nc[t₁]
	μ_new_f, d_new_f[f, t₁]
	μ_used_f, p_used_f[f, t₁]
	μ_new, d_new[t₁]
	μ_used, d_used[t₁]
end

# ==============================================================================
# Solve
# ==============================================================================
baseline = solve(calibration, db; replace_nothing=1.0)

println("Calibrated parameters:")
println("  μ_d = ", baseline[μ_d])
println("  μ_nc = ", baseline[μ_nc])
println("  μ_new = ", baseline[μ_new])
println("  μ_used = ", baseline[μ_used])
println("  μ_new_f = ", [baseline[μ_new_f[ff]] for ff in f])
println("  μ_used_f = ", [baseline[μ_used_f[ff]] for ff in f])

# ==============================================================================
# Counterfactual: EV subsidy (lower electric car price)
# ==============================================================================
T = 2040
t₁ = 2026

scenario_equations = @block db begin
	d[t = t₁:T],
	d[t] == μ_d * C[t] * (p_uc_d[t] / p_C[t])^(-σ_C)

	c_nc[t = t₁:T],
	c_nc[t] == μ_nc * C[t] * (p_nc[t] / p_C[t])^(-σ_C)

	p_C[t = t₁:T],
	p_C[t] * C[t] == p_uc_d[t] * d[t] + p_nc[t] * c_nc[t]

	d_new[t = t₁:T],
	d_new[t] == μ_new * d[t] * (p_uc_new[t] / p_d[t])^(-σ)

	d_used[t = t₁:T],
	d_used[t] == μ_used * d[t] * (p_uc_used[t] / p_d[t])^(-σ)

	p_d[t = t₁:T],
	p_d[t] * d[t] == p_uc_new[t] * d_new[t] + p_uc_used[t] * d_used[t]

	p_uc_d[t = t₁:T-1],
	p_uc_d[t] == p_d[t] - (1 - δ_d[t+1]) / (1 + r[t+1]) * p_d[t+1]

	p_uc_d[t = [T]],
	p_uc_d[t] == p_d[t] * (1 - (1 - δ_d[t]) / (1 + r[t]))

	d_new_f[f = f, t = t₁:T],
	d_new_f[f,t] == μ_new_f[f] * d_new[t] * (p_uc_new_f[f,t] / p_uc_new[t])^(-σ_new)

	p_uc_new[t = t₁:T],
	p_uc_new[t] * d_new[t] == ∑(p_uc_new_f[f,t] * d_new_f[f,t] for f ∈ f)

	p_uc_new_f[f = f, t = t₁:T-1],
	p_uc_new_f[f,t] == p_new_f[f,t] - (1 - δ0_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1]

	p_uc_new_f[f = f, t = [T]],
	p_uc_new_f[f,t] == p_new_f[f,t] - (1 - δ0_f[f,t]) / (1 + r[t]) * p_used_f[f,t]

	d_used_f[f = f, t = t₁:T],
	d_used_f[f,t] == (1 - δ_f[f,t]) * d_used_f[f,t-1] + (1 - δ0_f[f,t]) * d_new_f[f,t-1]

	p_used_f[f = f, t = t₁:T],
	d_used_f[f,t] - h_f[f] * d_used_f[f,t-1] == μ_used_f[f] * d_used[t] * (p_uc_f[f,t] / p_uc_used[t])^(-σ_used)

	p_uc_used[t = t₁:T],
	p_uc_used[t] * d_used[t] == ∑(p_uc_f[f,t] * (d_used_f[f,t] - h_f[f] * d_used_f[f,t-1]) for f ∈ f)

	p_uc_f[f = f, t = t₁:T-1],
	p_uc_f[f,t] == p_used_f[f,t] - (1 - δ_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1] + (1 - δ_f[f,t+1]) / (1 + r[t+1]) * h_f[f] * p_uc_used[t+1] * μ_used_f[f] * (d_used[t+1] / (d_used_f[f,t+1] - h_f[f] * d_used_f[f,t]))^(1/σ_used)

	p_uc_f[f = f, t = [T]],
	p_uc_f[f,t] == p_used_f[f,t] - (1 - δ_f[f,t]) / (1 + r[t]) * p_used_f[f,t] + (1 - δ_f[f,t]) / (1 + r[t]) * h_f[f] * p_uc_used[t] * μ_used_f[f] * (d_used[t] / (d_used_f[f,t] - h_f[f] * d_used_f[f,t-1]))^(1/σ_used)
end

scenario = copy(baseline)

for tt in t₁:T
	scenario[p_new_f[:electric, tt]] = 0.8  # 20% price reduction
end

solve!(scenario_equations, scenario)

multipliers = scenario ./ baseline .- 1

# ==============================================================================
# Plotting
# ==============================================================================
years = t₁:T

get_ts(data, var, yrs=years) = [data[var[tt]] for tt in yrs]
get_ts(data, var, ff::Symbol, yrs=years) = [data[var[ff, tt]] for tt in yrs]

pct(x) = x .* 100  # convert multiplier (fraction) to percentage points

fig = Figure(size=(1200, 900))

# --- New cars by fuel type ---
ax1 = Axis(fig[1,1]; title="New cars by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax1, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax1, years, pct(get_ts(multipliers, d_new_f, ff)); linestyle=ls, label=string(ff))
end
axislegend(ax1; position=:rt)

# --- Used-car stock by fuel type ---
ax2 = Axis(fig[1,2]; title="Used-car stock by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax2, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax2, years, pct(get_ts(multipliers, d_used_f, ff)); linestyle=ls, label=string(ff))
end
axislegend(ax2; position=:rt)

# --- Aggregates: new, used, total cars, non-car ---
ax3 = Axis(fig[2,1]; title="Aggregates", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax3, 0; color=:black, linewidth=0.8)
for (var, name_str, ls) in [(d_new, "New", :solid), (d_used, "Used", :dash), (d, "Cars total", :dot), (c_nc, "Non-car", :dashdot)]
	lines!(ax3, years, pct(get_ts(multipliers, var)); linestyle=ls, label=name_str)
end
axislegend(ax3; position=:rt)

# --- Prices: user costs ---
ax4 = Axis(fig[2,2]; title="User costs by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax4, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(get_ts(multipliers, p_uc_new_f, ff)); linestyle=ls, label="New $(ff)")
end
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(get_ts(multipliers, p_uc_f, ff)); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax4; position=:rt, nbanks=2)

Label(fig[0, :], "EV Subsidy Scenario: 20% Price Reduction (2026–$T)"; fontsize=20, font=:bold)

save("cars_scenario.svg", fig)
println("\nPlot saved to cars_scenario.svg")
