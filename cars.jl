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
using JuMP: Model, set_silent, name
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
const max_T = 2100
const plot_T = 2060
const t = t₀:max_T

t₁::Int = t₀ + 1
T::Int = max_T

# ==============================================================================
# Variables
# ==============================================================================
@variables db.model begin
	# Quantities — top level
	C[t], "Total consumption (exogenous)"
	C_noncar[t], "Non-car consumption"
	car[t], "Car-service aggregate"
	car_new[t], "New-car aggregate"
	car_used[t], "Used-car aggregate"
	car_new_f[f,t], "New cars by fuel type"
	car_used_f[f,t], "Used-car stock by fuel type"

	# Prices
	p_C[t], "Price of total consumption"
	p_nc[t], "Price of non-car consumption (exogenous)"
	p_car[t], "User cost of car-service aggregate"
	p_uc_new[t], "User cost of new-car aggregate"
	p_uc_used[t], "User cost of used-car aggregate"
	p_new_f[f,t], "Purchase price of new car by fuel type (exogenous)"
	p_used_f[f,t], "Price of used car by fuel type"
	p_uc_new_f[f,t], "User cost of new car by fuel type"
	p_uc_used_f[f,t], "User cost of used car by fuel type"

	# Exogenous
	r[t], "Interest rate"
	δ_new_f[f,t], "First-period depreciation (new→used)"
	δ_used_f[f,t], "Ongoing depreciation of used cars"

	# Parameters / calibrated
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

	# Tax/subsidy instruments
	τ_f[f,t], "Ad-valorem tax (+) or subsidy (-) on new cars by fuel type"
end

# ==============================================================================
# Equations
# ==============================================================================
equations() = @block db begin
	# ------------------------------------------------------------------
	# Top level: car aggregate vs non-car consumption
	# ------------------------------------------------------------------
	car[t = t₁:T],
	car[t] == μ_car * C[t] * (p_car[t] / p_C[t])^(-σ_C)

	C_noncar[t = t₁:T],
	C_noncar[t] == μ_noncar * C[t] * (p_nc[t] / p_C[t])^(-σ_C)

	p_C[t = t₁:T],
	p_C[t] * C[t] == p_car[t] * car[t] + p_nc[t] * C_noncar[t]

	# ------------------------------------------------------------------
	# Car nest: new vs used
	# ------------------------------------------------------------------
	car_new[t = t₁:T],
	car_new[t] == μ_new * car[t] * (p_uc_new[t] / p_car[t])^(-σ)

	car_used[t = t₁:T],
	car_used[t] == μ_used * car[t] * (p_uc_used[t] / p_car[t])^(-σ)

	p_car[t = t₁:T],
	p_car[t] * car[t] == p_uc_new[t] * car_new[t] + p_uc_used[t] * car_used[t]

	# ------------------------------------------------------------------
	# New-car nest: across fuel types
	# ------------------------------------------------------------------
	car_new_f[f = f, t = t₁:T],
	car_new_f[f,t] == μ_new_f[f] * car_new[t] * (p_uc_new_f[f,t] / p_uc_new[t])^(-σ_new)

	p_uc_new[t = t₁:T],
	p_uc_new[t] * car_new[t] == ∑(p_uc_new_f[f,t] * car_new_f[f,t] for f ∈ f)

	p_uc_new_f[f = f, t = t₁:T-1],
	p_uc_new_f[f,t] == p_new_f[f,t] * (1 + τ_f[f,t]) - (1 - δ_new_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1]

	p_uc_new_f[f = f, t = [T]],
	p_uc_new_f[f,t] == p_new_f[f,t] * (1 + τ_f[f,t]) - (1 - δ_new_f[f,t]) / (1 + r[t]) * p_used_f[f,t]

	# ------------------------------------------------------------------
	# Used-car nest: across fuel types (with habits)
	# ------------------------------------------------------------------
	car_used_f[f = f, t = t₁:T],
	car_used_f[f,t] == (1 - δ_used_f[f,t]) * car_used_f[f,t-1] + (1 - δ_new_f[f,t]) * car_new_f[f,t-1]

	# Price adjusts to clear the used-car market given stocks
	p_used_f[f = f, t = t₁:T],
	car_used_f[f,t] - h_f[f] * car_used_f[f,t-1] == μ_used_f[f] * car_used[t] * (p_uc_used_f[f,t] / p_uc_used[t])^(-σ_used)

	p_uc_used[t = t₁:T],
	p_uc_used[t] * car_used[t] == ∑(p_uc_used_f[f,t] * (car_used_f[f,t] - h_f[f] * car_used_f[f,t-1]) for f ∈ f)

	p_uc_used_f[f = f, t = t₁:T-1],
	p_uc_used_f[f,t] == p_used_f[f,t] - (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * p_used_f[f,t+1] + (1 - δ_used_f[f,t+1]) / (1 + r[t+1]) * h_f[f] * p_uc_used[t+1] * μ_used_f[f] * (car_used[t+1] / (car_used_f[f,t+1] - h_f[f] * car_used_f[f,t]))^(1/σ_used)

	p_uc_used_f[f = f, t = [T]],
	p_uc_used_f[f,t] == p_used_f[f,t] - (1 - δ_used_f[f,t]) / (1 + r[t]) * p_used_f[f,t] + (1 - δ_used_f[f,t]) / (1 + r[t]) * h_f[f] * p_uc_used[t] * μ_used_f[f] * (car_used[t] / (car_used_f[f,t] - h_f[f] * car_used_f[f,t-1]))^(1/σ_used)
end

# ==============================================================================
# Data
# ==============================================================================
# Elasticities
db[σ_C] = 0.5
db[σ] = 1.5
db[σ_new] = 3.0
db[σ_used] = 3.0

# Habits:
db[h_f] = [0.8, 0.8]

# Interest rate and depreciation
db[r] .= 0.04
db[δ_new_f] .= 0.25
db[δ_used_f] .= 0.10

# Exogenous
db[C] .= 1.0
db[p_nc] .= 1.0
db[p_new_f] .= 1.0
db[τ_f] .= 0.0

# ==============================================================================
# Calibration targets (equal shares)
# ==============================================================================
db[car_new_f[:, t₁]] .= 0.5 / length(f)
db[car_new[t₁]] = 0.5
db[car_used[t₁]] = 0.5
db[car[t₁]] = 0.5
db[C_noncar[t₁]] = 0.5
db[p_used_f[:, t₁]] .= 0.5

# Steady-state: lagged values equal t₁ values
steady_state() = @block db begin
	car_new_f[f = f, t = [t₀]],
	car_new_f[f,t] == car_new_f[f,t₁]

	car_used_f[f = f, t = [t₀]],
	car_used_f[f,t] == car_used_f[f,t₁]
end

# ==============================================================================
# Solve: calibrate steady state (T = t₁), then solve full horizon
# ==============================================================================
T = t₁
calibration = equations() + steady_state()
@endo_exo_swap! calibration begin
	μ_car, car[t₁]
	μ_noncar, C_noncar[t₁]
	μ_new_f, car_new_f[f, t₁]
	μ_used_f, p_used_f[f, t₁]
	μ_new, car_new[t₁]
	μ_used, car_used[t₁]
end
baseline = solve(calibration, db; replace_nothing=1.0)

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

T = max_T
extend_to_horizon!(baseline, calibration, t₁, t₁+1:T)
full_model = equations()
baseline = solve(full_model, baseline)

println("Calibrated parameters:")
println("  μ_car = ", baseline[μ_car])
println("  μ_noncar = ", baseline[μ_noncar])
println("  μ_new = ", baseline[μ_new])
println("  μ_used = ", baseline[μ_used])
println("  μ_new_f = ", [baseline[μ_new_f[ff]] for ff in f])
println("  μ_used_f = ", [baseline[μ_used_f[ff]] for ff in f])

# ==============================================================================
# Baseline sanity-check plots
# ==============================================================================
base_years = t₁:plot_T

fig_base = Figure(size=(1200, 900))

ax1 = Axis(fig_base[1,1]; title="Quantities by fuel type", xlabel="Year")
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax1, base_years, Float64.(baseline[car_new_f[ff, base_years]]); linestyle=ls, label="New $(ff)")
	lines!(ax1, base_years, Float64.(baseline[car_used_f[ff, base_years]]); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig_base[1,2]; title="Aggregates", xlabel="Year")
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax2, base_years, Float64.(baseline[var[base_years]]); linestyle=ls, label=name_str)
end
axislegend(ax2; position=:rt)

ax3 = Axis(fig_base[2,1]; title="Prices (user costs)", xlabel="Year")
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax3, base_years, Float64.(baseline[p_uc_new_f[ff, base_years]]); linestyle=ls, label="New $(ff)")
	lines!(ax3, base_years, Float64.(baseline[p_uc_used_f[ff, base_years]]); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
lines!(ax3, base_years, Float64.(baseline[p_car[base_years]]); linestyle=:dot, color=:black, label="p_car")
axislegend(ax3; position=:rt)

ax4 = Axis(fig_base[2,2]; title="Budget identity check", xlabel="Year")
budget_lhs = Float64.(baseline[p_C[base_years]]) .* Float64.(baseline[C[base_years]])
budget_rhs = Float64.(baseline[p_car[base_years]]) .* Float64.(baseline[car[base_years]]) .+ Float64.(baseline[p_nc[base_years]]) .* Float64.(baseline[C_noncar[base_years]])
lines!(ax4, base_years, budget_lhs; label="p_C · C")
lines!(ax4, base_years, budget_rhs; linestyle=:dash, label="p_car · car + p_nc · C_nc")
axislegend(ax4; position=:rt)

Label(fig_base[0, :], "Baseline Sanity Checks"; fontsize=20, font=:bold)

save("cars_baseline.svg", fig_base)
println("\nBaseline plot saved to cars_baseline.svg")

# ==============================================================================
# Counterfactual: EV subsidy (lower electric car price)
# ==============================================================================
T = max_T
t₁ = 2026
scenario = copy(baseline)
scenario[p_new_f[:electric, t₁:T]] .= 0.9  # 10% price reduction
solve!(equations(), scenario)

multipliers = scenario ./ baseline .- 1

# ==============================================================================
# Plotting
# ==============================================================================
years = t₁-1:plot_T

pct(x) = x .* 100  # convert multiplier (fraction) to percentage points

fig = Figure(size=(1200, 1200))

# --- New cars by fuel type ---
ax1 = Axis(fig[1,1]; title="New cars by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax1, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax1, years, pct(Float64.(multipliers[car_new_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax1; position=:rt)

# --- Used-car stock by fuel type ---
ax2 = Axis(fig[1,2]; title="Used-car stock by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax2, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax2, years, pct(Float64.(multipliers[car_used_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax2; position=:rt)

# --- Aggregates: new, used, total cars, non-car ---
ax3 = Axis(fig[2,1]; title="Aggregates", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax3, 0; color=:black, linewidth=0.8)
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars total", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax3, years, pct(Float64.(multipliers[var[years]])); linestyle=ls, label=name_str)
end
axislegend(ax3; position=:rt)

# --- Prices: user costs ---
ax4 = Axis(fig[2,2]; title="User costs by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax4, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers[p_uc_new_f[ff, years]])); linestyle=ls, label="New $(ff)")
end
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers[p_uc_used_f[ff, years]])); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax4; position=:rt, nbanks=2)

# --- Prices: used-car spot prices ---
ax5 = Axis(fig[3,1]; title="Used-car spot prices by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax5, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax5, years, pct(Float64.(multipliers[p_used_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax5; position=:rt)

# --- Prices: new-car purchase prices (exogenous shock) ---
ax6 = Axis(fig[3,2]; title="New-car purchase prices by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax6, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax6, years, pct(Float64.(multipliers[p_new_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax6; position=:rt)

Label(fig[0, :], "EV Subsidy Scenario: 20% Price Reduction ($(t₁)–$(T))"; fontsize=20, font=:bold)

save("cars_scenario1.svg", fig)
println("\nPlot saved to cars_scenario1.svg")

# ==============================================================================
# Counterfactual 2: Revenue-neutral tax/subsidy
# A 10% subsidy on electric cars financed by a constant endogenous tax on
# petrol cars, such that the present value of net tax revenue is zero.
# ==============================================================================
@variables db.model begin
	tax_revenue[t], "Net tax revenue from new-car taxes/subsidies"
	pv_tax_revenue, "Present value of net tax revenue"
end

revenue_neutral_eqs() = @block db begin
	tax_revenue[t = t₁:T],
	tax_revenue[t] == ∑(τ_f[f,t] * p_new_f[f,t] * car_new_f[f,t] for f ∈ f)

	pv_tax_revenue,
	pv_tax_revenue == ∑(tax_revenue[t] / (1 + r[t])^(t - t₁) for t ∈ t₁:T)

	τ_f[f = [:petrol], t = t₁:T-1],
	τ_f[f, t] == τ_f[f, T]
end

scenario2 = copy(baseline)
scenario2[τ_f[:electric, t₁:T]] .= -0.10  # 10% subsidy on electric
scenario2[τ_f[:petrol, t₁:T]] .= 0.3      # initial guess
scenario2[tax_revenue[t₁:T]] .= 0.0
scenario2[pv_tax_revenue] = 0.0

model2 = equations() + revenue_neutral_eqs()
@endo_exo_swap! model2 begin
	τ_f[:petrol, T], pv_tax_revenue
end
solve!(model2, scenario2)

τ_val = round(scenario2[τ_f[:petrol, t₁]] * 100; digits=2)
println("\nRevenue-neutral petrol tax rate: τ_petrol = $τ_val%")

# ==============================================================================
# Plotting — Revenue-neutral scenario
# ==============================================================================
add_missing_model_variables!(baseline)
multipliers2 = scenario2 ./ baseline .- 1

fig2 = Figure(size=(1200, 1200))

ax1 = Axis(fig2[1,1]; title="New cars by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax1, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax1, years, pct(Float64.(multipliers2[car_new_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig2[1,2]; title="Used-car stock by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax2, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax2, years, pct(Float64.(multipliers2[car_used_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax2; position=:rt)

ax3 = Axis(fig2[2,1]; title="Aggregates", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax3, 0; color=:black, linewidth=0.8)
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars total", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax3, years, pct(Float64.(multipliers2[var[years]])); linestyle=ls, label=name_str)
end
axislegend(ax3; position=:rt)

ax4 = Axis(fig2[2,2]; title="User costs by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax4, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers2[p_uc_new_f[ff, years]])); linestyle=ls, label="New $(ff)")
end
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers2[p_uc_used_f[ff, years]])); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax4; position=:rt, nbanks=2)

ax5 = Axis(fig2[3,1]; title="Used-car spot prices by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax5, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax5, years, pct(Float64.(multipliers2[p_used_f[ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax5; position=:rt)

ax6 = Axis(fig2[3,2]; title="Implied ad-valorem tax/subsidy rates", xlabel="Year", ylabel="Rate")
hlines!(ax6, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax6, years, Float64.(scenario2[τ_f[ff, years]]); linestyle=ls, label=string(ff))
end
axislegend(ax6; position=:rt)

Label(fig2[0, :], "PV-Neutral Tax/Subsidy: 10% EV Subsidy ($(t₁)–$(T))"; fontsize=20, font=:bold)

save("cars_scenario2.svg", fig2)
println("\nPlot saved to cars_scenario2.svg")
