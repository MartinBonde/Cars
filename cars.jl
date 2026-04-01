# Partial Equilibrium Model of New and Used Cars — Scenario Analysis
#
# Calibrates the model from car_model.jl, then runs two counterfactuals:
#   1. EV subsidy (10% reduction in electric car purchase price)
#   2. Revenue-neutral tax/subsidy (10% EV subsidy financed by petrol tax)

include("car_model.jl")

const plot_T = 2060

# ==============================================================================
# Calibrate baseline (σ = 1.5 for new-vs-used; β_h = 1.0)
# ==============================================================================
baseline = calibrate()

# ==============================================================================
# Baseline sanity-check plots
# ==============================================================================
base_years = t₁:plot_T

fig_base = Figure(size=(1200, 900))

# Aggregate new/used by fuel type for plotting: sum over brands
new_by_f(data, ff, yrs) = Float64.([sum(data[car_new_bf[bb, ff, yr]] for bb in b) for yr in yrs])
used_by_f(data, ff, yrs) = Float64.([sum(data[car_used_bf[bb, ff, yr]] for bb in b) for yr in yrs])

ax1 = Axis(fig_base[1,1]; title="Quantities by fuel type", xlabel="Year")
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax1, base_years, new_by_f(baseline, ff, base_years); linestyle=ls, label="New $(ff)")
	lines!(ax1, base_years, used_by_f(baseline, ff, base_years); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig_base[1,2]; title="Aggregates", xlabel="Year")
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax2, base_years, Float64.(baseline[var[base_years]]); linestyle=ls, label=name_str)
end
axislegend(ax2; position=:rt)

# Use brand1 user costs as representative (symmetric brands)
ax3 = Axis(fig_base[2,1]; title="Prices (user costs, brand1)", xlabel="Year")
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax3, base_years, Float64.(baseline[p_uc_new_bf[:brand1, ff, base_years]]); linestyle=ls, label="New $(ff)")
	lines!(ax3, base_years, Float64.(baseline[p_uc_used_bf[:brand1, ff, base_years]]); linestyle=ls, color=Cycled(3), label="Used $(ff)")
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
# Counterfactual 1: EV subsidy (lower electric car price for all brands)
# ==============================================================================
T = max_T
t₁ = 2026
scenario = copy(baseline)
scenario[p_new_bf[:, :electric, t₁:T]] .= 0.9  # 10% price reduction
solve!(equations(), scenario)

multipliers = scenario ./ baseline .- 1

# ==============================================================================
# Plotting — EV subsidy
# ==============================================================================
years = t₁-1:plot_T

pct(x) = x .* 100

fig = Figure(size=(1200, 1200))

ax1 = Axis(fig[1,1]; title="New cars by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax1, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	new_base = new_by_f(baseline, ff, years)
	new_scen = new_by_f(scenario, ff, years)
	lines!(ax1, years, pct(new_scen ./ new_base .- 1); linestyle=ls, label=string(ff))
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig[1,2]; title="Used-car stock by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax2, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	used_base = used_by_f(baseline, ff, years)
	used_scen = used_by_f(scenario, ff, years)
	lines!(ax2, years, pct(used_scen ./ used_base .- 1); linestyle=ls, label=string(ff))
end
axislegend(ax2; position=:rt)

ax3 = Axis(fig[2,1]; title="Aggregates", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax3, 0; color=:black, linewidth=0.8)
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars total", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax3, years, pct(Float64.(multipliers[var[years]])); linestyle=ls, label=name_str)
end
axislegend(ax3; position=:rt)

ax4 = Axis(fig[2,2]; title="User costs by fuel type (brand1)", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax4, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers[p_uc_new_bf[:brand1, ff, years]])); linestyle=ls, label="New $(ff)")
end
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers[p_uc_used_bf[:brand1, ff, years]])); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax4; position=:rt, nbanks=2)

ax5 = Axis(fig[3,1]; title="Used-car spot prices by fuel type (brand1)", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax5, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax5, years, pct(Float64.(multipliers[p_used_bf[:brand1, ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax5; position=:rt)

ax6 = Axis(fig[3,2]; title="New-car purchase prices by fuel type (brand1)", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax6, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax6, years, pct(Float64.(multipliers[p_new_bf[:brand1, ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax6; position=:rt)

Label(fig[0, :], "EV Subsidy Scenario: 10% Price Reduction ($(t₁)–$(T))"; fontsize=20, font=:bold)

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
	tax_revenue[t] == ∑(τ_bf[b,f,t] * p_new_bf[b,f,t] * car_new_bf[b,f,t] for b ∈ b, f ∈ f)

	pv_tax_revenue,
	pv_tax_revenue == ∑(tax_revenue[t] / (1 + r[t])^(t - t₁) for t ∈ t₁:T)

	τ_bf[b = b, f = [:petrol], t = t₁:T-1],
	τ_bf[b, f, t] == τ_bf[b, f, T]
end

scenario2 = copy(baseline)
scenario2[τ_bf[:, :electric, t₁:T]] .= -0.10
scenario2[τ_bf[:, :petrol, t₁:T]] .= 0.3
scenario2[tax_revenue[t₁:T]] .= 0.0
scenario2[pv_tax_revenue] = 0.0

model2 = equations() + revenue_neutral_eqs()
@endo_exo_swap! model2 begin
	τ_bf[:brand1, :petrol, T], pv_tax_revenue
end
solve!(model2, scenario2)

τ_val = round(scenario2[τ_bf[:brand1, :petrol, t₁]] * 100; digits=2)
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
	new_base = new_by_f(baseline, ff, years)
	new_scen = new_by_f(scenario2, ff, years)
	lines!(ax1, years, pct(new_scen ./ new_base .- 1); linestyle=ls, label=string(ff))
end
axislegend(ax1; position=:rt)

ax2 = Axis(fig2[1,2]; title="Used-car stock by fuel type", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax2, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	used_base = used_by_f(baseline, ff, years)
	used_scen = used_by_f(scenario2, ff, years)
	lines!(ax2, years, pct(used_scen ./ used_base .- 1); linestyle=ls, label=string(ff))
end
axislegend(ax2; position=:rt)

ax3 = Axis(fig2[2,1]; title="Aggregates", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax3, 0; color=:black, linewidth=0.8)
for (var, name_str, ls) in [(car_new, "New", :solid), (car_used, "Used", :dash), (car, "Cars total", :dot), (C_noncar, "Non-car", :dashdot)]
	lines!(ax3, years, pct(Float64.(multipliers2[var[years]])); linestyle=ls, label=name_str)
end
axislegend(ax3; position=:rt)

ax4 = Axis(fig2[2,2]; title="User costs by fuel type (brand1)", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax4, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers2[p_uc_new_bf[:brand1, ff, years]])); linestyle=ls, label="New $(ff)")
end
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax4, years, pct(Float64.(multipliers2[p_uc_used_bf[:brand1, ff, years]])); linestyle=ls, color=Cycled(3), label="Used $(ff)")
end
axislegend(ax4; position=:rt, nbanks=2)

ax5 = Axis(fig2[3,1]; title="Used-car spot prices by fuel type (brand1)", xlabel="Year", ylabel="% deviation from baseline")
hlines!(ax5, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax5, years, pct(Float64.(multipliers2[p_used_bf[:brand1, ff, years]])); linestyle=ls, label=string(ff))
end
axislegend(ax5; position=:rt)

ax6 = Axis(fig2[3,2]; title="Implied ad-valorem tax/subsidy rates", xlabel="Year", ylabel="Rate")
hlines!(ax6, 0; color=:black, linewidth=0.8)
for (ff, ls) in zip(f, [:solid, :dash])
	lines!(ax6, years, Float64.(scenario2[τ_bf[:brand1, ff, years]]); linestyle=ls, label=string(ff))
end
axislegend(ax6; position=:rt)

Label(fig2[0, :], "PV-Neutral Tax/Subsidy: 10% EV Subsidy ($(t₁)–$(T))"; fontsize=20, font=:bold)

save("cars_scenario2.svg", fig2)
println("\nPlot saved to cars_scenario2.svg")
