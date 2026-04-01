# Effect of the Number of Brands on Market Power
#
# For each brand count N, spawns a separate Julia process with a patched
# car_model.jl (different const b), calibrates, and computes brand-level
# demand elasticities at β_h = 1 and β_h = 0.

using CairoMakie

const brand_counts = [1, 2, 3, 5, 8, 12]

worker_script = raw"""
# Patched car_model.jl is included via BRAND_MODEL_PATH
include(ENV["BRAND_MODEL_PATH"])

const ε = 0.01
const short_T = 2040
const shock_brand = first(b)

baseline_ss = calibrate(full_horizon=false)

function _ss_elast(baseline_ss, β_h_val)
	cf = copy(baseline_ss)
	cf[β_h] = β_h_val
	global T = t₁
	ss_model = equations() + steady_state()
	solve!(ss_model, cf)
	shocked = copy(cf)
	shocked[p_new_bf[shock_brand, :, t₁]] .= cf[p_new_bf[shock_brand, :, t₁]] .* (1 + ε)
	solve!(ss_model, shocked)
	global T = max_T
	return (shocked[car_new_b[shock_brand, t₁]] / cf[car_new_b[shock_brand, t₁]] - 1) * 100
end

function _imp_elast(baseline_ss, β_h_val)
	cf = copy(baseline_ss)
	cf[β_h] = β_h_val
	global T = t₁
	ss_model = equations() + steady_state()
	solve!(ss_model, cf)
	global T = short_T
	extend_to_horizon!(cf, ss_model, t₁, t₁+1:T)
	fm = equations()
	solve!(fm, cf)
	shocked = copy(cf)
	shock_year = 2026
	shocked[p_new_bf[shock_brand, :, shock_year:T]] .= cf[p_new_bf[shock_brand, :, shock_year:T]] .* (1 + ε)
	solve!(fm, shocked)
	q = (shocked[car_new_b[shock_brand, shock_year]] / cf[car_new_b[shock_brand, shock_year]] - 1) * 100
	global T = max_T
	return q
end

imp_fwd = _imp_elast(baseline_ss, 1.0)
ss_fwd = _ss_elast(baseline_ss, 1.0)
imp_myopic = _imp_elast(baseline_ss, 0.0)
ss_myopic = _ss_elast(baseline_ss, 0.0)

println("RESULT:$(imp_fwd):$(ss_fwd):$(imp_myopic):$(ss_myopic)")
"""

results = Dict{Int, NamedTuple{(:imp_fwd, :ss_fwd, :imp_myopic, :ss_myopic), NTuple{4, Float64}}}()
base_dir = @__DIR__

for N in brand_counts
	println("\n" * "=" ^ 60)
	println("N = $N brands")
	println("=" ^ 60)

	brand_syms = [":brand$i" for i in 1:N]
	brand_line = "const b = [" * join(brand_syms, ", ") * "]"

	car_model_src = read(joinpath(base_dir, "car_model.jl"), String)
	patched_src = replace(car_model_src, r"const b = \[[^\n]*\]" => brand_line)

	patched_model_path = joinpath(tempdir(), "car_model_N$(N).jl")
	write(patched_model_path, patched_src)

	worker_path = joinpath(tempdir(), "worker_N$(N).jl")
	write(worker_path, worker_script)

	cmd = setenv(`julia --project=$(base_dir) $(worker_path)`,
		"BRAND_MODEL_PATH" => patched_model_path)
	output = read(cmd, String)

	for line in split(output, '\n')
		if startswith(line, "RESULT:")
			parts = split(line, ':')
			r = (imp_fwd=parse(Float64, parts[2]),
			     ss_fwd=parse(Float64, parts[3]),
			     imp_myopic=parse(Float64, parts[4]),
			     ss_myopic=parse(Float64, parts[5]))
			results[N] = r
			println("  fwd:    impact=$(round(r.imp_fwd, digits=4))%  ss=$(round(r.ss_fwd, digits=4))%")
			println("  myopic: impact=$(round(r.imp_myopic, digits=4))%  ss=$(round(r.ss_myopic, digits=4))%")
		end
	end
end

# ==============================================================================
# Print summary
# ==============================================================================
println("\n" * "=" ^ 60)
println("Summary: Effect of Number of Brands")
println("=" ^ 60)
Ns = sort(collect(keys(results)))
for N in Ns
	r = results[N]
	gap = r.imp_myopic - r.imp_fwd
	println("  N=$(lpad(N,2)):  imp_fwd=$(round(r.imp_fwd, digits=3))  imp_myopic=$(round(r.imp_myopic, digits=3))  gap=$(round(gap, digits=3))  ss_fwd=$(round(r.ss_fwd, digits=3))  ss_myopic=$(round(r.ss_myopic, digits=3))")
end

# ==============================================================================
# Plot
# ==============================================================================
imp_fwd = [results[n].imp_fwd for n in Ns]
imp_myopic = [results[n].imp_myopic for n in Ns]

c_fwd    = RGBAf(0.122, 0.467, 0.706, 1.0)
c_myopic = RGBAf(0.890, 0.467, 0.157, 1.0)
c_gap    = RGBAf(0.5, 0.2, 0.6, 0.12)
c_base   = RGBAf(0.15, 0.15, 0.15, 1.0)

fig = Figure(size=(800, 500), backgroundcolor=:white)
ax = Axis(fig[1,1];
	title="Brand-Level Market Power vs Number of Brands  (σ_brand = 5)",
	xlabel="Number of symmetric brands",
	ylabel="Δ brand quantity  (% per 1 % cost-push)",
	titlefont=:bold, titlesize=16)
ax.xgridvisible = false
ax.ygridvisible = true
ax.ygridcolor = RGBAf(0, 0, 0, 0.06)
ax.spinewidth = 0.8
ax.topspinevisible = false
ax.rightspinevisible = false

band!(ax, Float64.(Ns), imp_fwd, imp_myopic; color=c_gap)
lines!(ax, Float64.(Ns), imp_fwd; color=c_fwd, linewidth=2.5, label="βₕ = 1 (forward-looking)")
scatter!(ax, Float64.(Ns), imp_fwd; color=c_fwd, markersize=7)
lines!(ax, Float64.(Ns), imp_myopic; color=c_myopic, linewidth=2.5, linestyle=:dash, label="βₕ = 0 (myopic)")
scatter!(ax, Float64.(Ns), imp_myopic; color=c_myopic, markersize=7)

bi = findfirst(==(5), Ns)
if bi !== nothing
	scatter!(ax, [5.0], [imp_fwd[bi]]; color=c_base, markersize=12, marker=:diamond)
	scatter!(ax, [5.0], [imp_myopic[bi]]; color=c_base, markersize=12, marker=:diamond)
end

axislegend(ax; position=:rb, framevisible=false, labelsize=11)

Label(fig[2,1],
	"Shaded area = forward-lookingness premium.  ◆ = baseline (5 brands).  Impact elasticity shown.";
	fontsize=10, halign=:center, color=:grey45, padding=(0, 0, 0, 2))
rowgap!(fig.layout, 1, 4)

save("brand_count.svg", fig)
println("\nPlot saved to brand_count.svg")
