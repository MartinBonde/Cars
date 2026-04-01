# Market Power of Individual Brand Owners
#
# Measures how a single brand owner's market power depends on:
#   a) Durability of used cars (ongoing depreciation δ)
#   b) Habit parameter h_f
#   c) Forward-lookingness wrt habits (β_h)
#   d) Brand substitutability (σ_brand) — does more variety strengthen
#      or weaken the forward-lookingness effect?
#
# The cost-push shock is applied to ONE brand only (brand1).
# The brand owner internalises within-brand fuel-type substitution but
# competes against other brands at the σ_brand level.

include("car_model.jl")

# ==============================================================================
# Calibrate ONCE at baseline parameters
# ==============================================================================
baseline_ss = calibrate(full_horizon=false)

# ==============================================================================
# Elasticity computation
# ==============================================================================
const ε = 0.01
const short_T = 2040
const shock_brand = :brand1

function is_valid(q_pct)
	abs(q_pct) < 10
end

"""Compute steady-state demand elasticity of a single brand to a brand-specific cost-push."""
function ss_elasticity(baseline_ss; δ_new=0.25, δ_used=0.10, h=0.8, β_h_val=1.0, σ_brand_val=5.0)
	cf = copy(baseline_ss)
	cf[δ_new_f] .= δ_new
	cf[δ_used_f] .= δ_used
	cf[h_f] = [h, h]
	cf[β_h] = β_h_val
	cf[σ_brand] = σ_brand_val

	global T = t₁
	ss_model = equations() + steady_state()
	solve!(ss_model, cf)

	shocked = copy(cf)
	shocked[p_new_bf[shock_brand, :, t₁]] .= cf[p_new_bf[shock_brand, :, t₁]] .* (1 + ε)
	solve!(ss_model, shocked)
	global T = max_T

	q_pct = (shocked[car_new_b[shock_brand, t₁]] / cf[car_new_b[shock_brand, t₁]] - 1) * 100
	is_valid(q_pct) || return nothing
	return q_pct
end

"""Compute impact (short-run dynamic) demand elasticity of a single brand."""
function impact_elasticity(baseline_ss; δ_new=0.25, δ_used=0.10, h=0.8, β_h_val=1.0, σ_brand_val=5.0, shock_year=2026)
	cf = copy(baseline_ss)
	cf[δ_new_f] .= δ_new
	cf[δ_used_f] .= δ_used
	cf[h_f] = [h, h]
	cf[β_h] = β_h_val
	cf[σ_brand] = σ_brand_val

	global T = t₁
	ss_model = equations() + steady_state()
	solve!(ss_model, cf)

	global T = short_T
	extend_to_horizon!(cf, ss_model, t₁, t₁+1:T)
	fm = equations()
	solve!(fm, cf)

	shocked = copy(cf)
	shocked[p_new_bf[shock_brand, :, shock_year:T]] .= cf[p_new_bf[shock_brand, :, shock_year:T]] .* (1 + ε)
	solve!(fm, shocked)
	global T = max_T

	q_pct = (shocked[car_new_b[shock_brand, shock_year]] / cf[car_new_b[shock_brand, shock_year]] - 1) * 100
	is_valid(q_pct) || return nothing
	return q_pct
end

# ==============================================================================
# a) Durability sweep (ongoing depreciation δ, δ₀ fixed)
# ==============================================================================
println("\n" * "=" ^ 60)
println("(a) Sweeping ongoing depreciation δ...")
println("=" ^ 60)

δ_values = [0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.25, 0.30]

dur_δ = Float64[]
dur_impact = Float64[]
dur_ss = Float64[]

for δ_u in δ_values
	print("  δ=$(round(δ_u, digits=3))... ")
	ss = ss_elasticity(baseline_ss; δ_used=δ_u)
	imp = impact_elasticity(baseline_ss; δ_used=δ_u)
	if ss !== nothing && imp !== nothing
		push!(dur_δ, δ_u)
		push!(dur_impact, imp)
		push!(dur_ss, ss)
		println("✓ impact=$(round(imp, digits=4))%, ss=$(round(ss, digits=4))%")
	else
		println("✗ solver anomaly, skipping")
	end
end

# ==============================================================================
# b) Habit sweep
# ==============================================================================
println("\n" * "=" ^ 60)
println("(b) Sweeping habit parameter h...")
println("=" ^ 60)

h_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.45, 0.55, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

hab_h = Float64[]
hab_impact = Float64[]
hab_ss = Float64[]

for hv in h_values
	print("  h=$(hv)... ")
	ss = ss_elasticity(baseline_ss; h=hv)
	imp = impact_elasticity(baseline_ss; h=hv)
	if ss !== nothing && imp !== nothing
		push!(hab_h, hv)
		push!(hab_impact, imp)
		push!(hab_ss, ss)
		println("✓ impact=$(round(imp, digits=4))%, ss=$(round(ss, digits=4))%")
	else
		println("✗ solver anomaly, skipping")
	end
end

# ==============================================================================
# c) Forward-lookingness sweep (β_h), h=0.8 fixed
# ==============================================================================
println("\n" * "=" ^ 60)
println("(c) Sweeping β_h (h=0.8 fixed)...")
println("=" ^ 60)

β_h_values = [0.0, 0.1, 0.15, 0.25, 0.3, 0.35, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

fwd_β = Float64[]
fwd_impact = Float64[]
fwd_ss = Float64[]

for bv in β_h_values
	print("  β_h=$(bv)... ")
	ss = ss_elasticity(baseline_ss; β_h_val=bv)
	imp = impact_elasticity(baseline_ss; β_h_val=bv)
	if ss !== nothing && imp !== nothing
		push!(fwd_β, bv)
		push!(fwd_impact, imp)
		push!(fwd_ss, ss)
		println("✓ impact=$(round(imp, digits=4))%, ss=$(round(ss, digits=4))%")
	else
		println("✗ solver anomaly, skipping")
	end
end

# ==============================================================================
# d) Brand substitutability × forward-lookingness interaction
#    For each σ_brand, compute the gap in brand-level elasticity between
#    β_h = 1 (fully forward-looking) and β_h = 0 (myopic).
# ==============================================================================
println("\n" * "=" ^ 60)
println("(d) Sweeping σ_brand (β_h = 0 vs β_h = 1)...")
println("=" ^ 60)

σ_brand_values = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0]

brand_σ = Float64[]
brand_impact_fwd = Float64[]
brand_impact_myopic = Float64[]
brand_ss_fwd = Float64[]
brand_ss_myopic = Float64[]

for σb in σ_brand_values
	print("  σ_brand=$(σb)... ")
	ss_f = ss_elasticity(baseline_ss; σ_brand_val=σb, β_h_val=1.0)
	ss_m = ss_elasticity(baseline_ss; σ_brand_val=σb, β_h_val=0.0)
	imp_f = impact_elasticity(baseline_ss; σ_brand_val=σb, β_h_val=1.0)
	imp_m = impact_elasticity(baseline_ss; σ_brand_val=σb, β_h_val=0.0)
	if all(x -> x !== nothing, (ss_f, ss_m, imp_f, imp_m))
		push!(brand_σ, σb)
		push!(brand_impact_fwd, imp_f)
		push!(brand_impact_myopic, imp_m)
		push!(brand_ss_fwd, ss_f)
		push!(brand_ss_myopic, ss_m)
		println("✓  fwd: imp=$(round(imp_f, digits=4))% ss=$(round(ss_f, digits=4))%  " *
		        "myopic: imp=$(round(imp_m, digits=4))% ss=$(round(ss_m, digits=4))%")
	else
		println("✗ solver anomaly, skipping")
	end
end

# ==============================================================================
# Plotting
# ==============================================================================
c_impact = RGBAf(0.122, 0.467, 0.706, 1.0)   # Tableau blue
c_ss     = RGBAf(0.890, 0.467, 0.157, 1.0)    # Tableau orange
c_fill   = RGBAf(0.122, 0.467, 0.706, 0.10)   # faint blue fill
c_base   = RGBAf(0.15, 0.15, 0.15, 1.0)       # near-black for baseline marker
c_fwd    = RGBAf(0.122, 0.467, 0.706, 1.0)    # blue for forward-looking
c_myopic = RGBAf(0.890, 0.467, 0.157, 1.0)    # orange for myopic
c_gap    = RGBAf(0.5, 0.2, 0.6, 0.12)         # faint purple fill for gap

function style_ax!(ax)
	ax.xgridvisible = false
	ax.ygridvisible = true
	ax.ygridcolor = RGBAf(0, 0, 0, 0.06)
	ax.spinewidth = 0.8
	ax.topspinevisible = false
	ax.rightspinevisible = false
	ax.xlabelsize = 13
	ax.ylabelsize = 13
	ax.titlesize = 14
	ax.titlefont = :bold
	ax.yticklabelsize = 11
	ax.xticklabelsize = 11
end

function plot_sweep!(ax, x, impact, ss; baseline_x=nothing)
	band!(ax, x, ss, impact; color=c_fill)
	lines!(ax, x, impact; color=c_impact, linewidth=2.5)
	scatter!(ax, x, impact; color=c_impact, markersize=5)
	lines!(ax, x, ss;     color=c_ss, linewidth=2.5, linestyle=:dash)
	scatter!(ax, x, ss;   color=c_ss, markersize=5)
	if baseline_x !== nothing
		bi = findfirst(==(baseline_x), x)
		if bi !== nothing
			scatter!(ax, [x[bi]], [impact[bi]]; color=c_base, markersize=12, marker=:diamond)
			scatter!(ax, [x[bi]], [ss[bi]];     color=c_base, markersize=12, marker=:diamond)
		end
	end
end

fig = Figure(size=(1500, 1000), backgroundcolor=:white)
ylabel_str = "Δ brand quantity  (% per 1 % cost-push)"

# ── Title row ────────────────────────────────────────────────────────────────
Label(fig[1, 1:2],
	"Market Power of Individual Brand Owners  ($(length(b)) symmetric brands)";
	fontsize=20, font=:bold, halign=:center, padding=(0, 0, 6, 0))

# ── Panel (a): Durability ────────────────────────────────────────────────────
ax1 = Axis(fig[2,1]; title="(a)  Durability",
	xlabel="Used-car survival rate  (1 − δ)", ylabel=ylabel_str)
style_ax!(ax1)
plot_sweep!(ax1, 1.0 .- dur_δ, dur_impact, dur_ss; baseline_x=0.90)

# ── Panel (b): Habit persistence ─────────────────────────────────────────────
ax2 = Axis(fig[2,2]; title="(b)  Habit persistence  (h)",
	xlabel="Habit parameter  h", ylabel="")
style_ax!(ax2)
ax2.ylabelvisible = false
plot_sweep!(ax2, hab_h, hab_impact, hab_ss; baseline_x=0.8)

# ── Panel (c): Forward-lookingness ───────────────────────────────────────────
ax3 = Axis(fig[3,1]; title="(c)  Forward-lookingness  (βₕ ,  h = 0.8)",
	xlabel="Habit-premium discount  βₕ\n(0 = myopic  ·  1 = fully forward-looking)", ylabel=ylabel_str)
style_ax!(ax3)
plot_sweep!(ax3, fwd_β, fwd_impact, fwd_ss; baseline_x=1.0)

# ── Panel (d): Brand variety × forward-lookingness ───────────────────────────
ax4 = Axis(fig[3,2]; title="(d)  Brand substitutability × forward-lookingness",
	xlabel="Brand elasticity  σ_brand\n(higher = more substitutable brands)", ylabel="")
style_ax!(ax4)
ax4.ylabelvisible = false

band!(ax4, brand_σ, brand_impact_fwd, brand_impact_myopic; color=c_gap)
lines!(ax4, brand_σ, brand_impact_fwd; color=c_fwd, linewidth=2.5, label="βₕ = 1 (forward-looking)")
scatter!(ax4, brand_σ, brand_impact_fwd; color=c_fwd, markersize=5)
lines!(ax4, brand_σ, brand_impact_myopic; color=c_myopic, linewidth=2.5, linestyle=:dash, label="βₕ = 0 (myopic)")
scatter!(ax4, brand_σ, brand_impact_myopic; color=c_myopic, markersize=5)
bi = findfirst(==(5.0), brand_σ)
if bi !== nothing
	scatter!(ax4, [brand_σ[bi]], [brand_impact_fwd[bi]]; color=c_base, markersize=12, marker=:diamond)
	scatter!(ax4, [brand_σ[bi]], [brand_impact_myopic[bi]]; color=c_base, markersize=12, marker=:diamond)
end
axislegend(ax4; position=:rb, framevisible=false, labelsize=10)

# ── Shared legend ────────────────────────────────────────────────────────────
Legend(fig[4, 1:2],
	[LineElement(color=c_impact, linewidth=2.5),
	 LineElement(color=c_ss, linewidth=2.5, linestyle=:dash),
	 [PolyElement(color=c_fill)],
	 MarkerElement(color=c_base, marker=:diamond, markersize=11)],
	["Impact (short run)", "Steady state (long run)",
	 "Short-run / long-run gap", "Baseline"];
	orientation=:horizontal, framevisible=false, labelsize=11,
	padding=(0, 0, 0, 0), colgap=24,
	tellwidth=false, tellheight=true, halign=:center)

Label(fig[5, 1:2],
	"More negative  ←  more elastic demand  =  less market power for brand owner.   " *
	"Share parameters calibrated once at baseline (◆).   " *
	"Shock applied to $(shock_brand) only.";
	fontsize=10, halign=:center, color=:grey45, padding=(0, 0, 0, 2))

# Tighten layout
rowgap!(fig.layout, 1, 2)
rowgap!(fig.layout, 2, 10)
rowgap!(fig.layout, 3, 6)
rowgap!(fig.layout, 4, 2)
colgap!(fig.layout, 16)

save("market_power.svg", fig)
println("\nPlot saved to market_power.svg")
