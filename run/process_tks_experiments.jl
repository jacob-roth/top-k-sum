import Pkg; Pkg.activate("."); Pkg.instantiate()
using LinearAlgebra, Random, DelimitedFiles, LaTeXStrings, LaTeXTabulars, Printf, NaNStatistics
global const PROJPATH = match(r".*top-k-sum/", @__DIR__).match
include(PROJPATH * "run/helper_tks_experiments.jl")
include(PROJPATH * "run/plot_tks_experiments.jl")
global const DATAPATH = "/home/roth0674/drive/tks_results/"
# global const DATAPATH = "/Users/jakeroth/Downloads/mks_results_20231007/"
global const DT = Float64
global const maxn_grid = 100_000
global const maxn_gurobi = 100_000

#
# inputs
#

out = load_tks_results(DATAPATH, 1, 100, maxn_grid, maxn_gurobi)
calc_bestfeasunsorted!(out)
nlevel =  10 .^ collect(1:7)
algl_names = ["ESGS", "PLCP", "GRID", "GRBS", "GRBU"]
algl_names_short = ["1", "2", "3", "4", "5"]
algl_names_detail = ["ESGS", "PLCP", "GRID", "GRBS", "GRBU"]
ri = [6,8,11,6,8,11]; out["rlevel"][ri]
ki = [2,2,2,4,4,4]; out["klevel"][ki]
nlevelj = [1_000, 10_000, 100_000, 1_000_000, 10_000_000]
nlevel_all = 10 .^ collect(1:7)
algl = [1,2,3,4,5]

out_large = load_tks_results(DATAPATH, 101, 102, 10^7, 10^7, true)

sci_1a(d) = begin
  """one-exponent-digit and one-base-decimaldigit"""
  expo = floor(log10(abs(d)))
  base = d * 10^-expo
  if d==0
    L"0"
  elseif isnan(base) || isinf(base) || isnan(expo) || isinf(expo)
    L"-"
  else
    @sprintf("%0.1fe%+d", base, expo)
  end
end
sci_1b(d) = begin
  """one-exponent-digit and zero-base-decimaldigit"""
  expo = floor(log10(abs(d)))
  base = d * 10^-expo
  if d==0
    L"0"
  elseif isnan(base) || isinf(base) || isnan(expo) || isinf(expo)
    L"-"
  else
    @sprintf("%0.0fe%d", base, expo)
  end
end
pct_safe(d) = begin
  if isnan(d) || isinf(d)
    L"-"
  else
    @sprintf("%0.2f", d)
  end
end

mnstd(x) = sci_1a(nanmean(x)) * " (" * sci_1b(nanstd(x)) * ")"
mn_std(x,y) = sci_1a(x) * " (" * sci_1b(y) * ")"

function sort_time_row(out::Dict, nlevelj::Vector)
  idx = out["sort_time"][:,1] .∈ Ref(nlevelj)
  data_mn = vec(nanmean(out["sort_time"][idx, 2:end], dims=2))
  ["Sort time (full)", [sci_1a(x) for x in data_mn]...]
end
function partial_sort_time_row(out::Dict, nlevelj::Vector)
  idx = out["partial_sort_time"][:,1] .∈ Ref(nlevelj)
  data_mn = vec(nanmean(out["partial_sort_time"][idx, 2:end], dims=2))
  ["Sort time (top-1\\%)", [sci_1a(x) for x in data_mn]...]
end

function experiment_time(exp_i::Integer,
  nlevelj::Vector,
  algl::Vector,
  algl_names::Vector,
  ri::Vector,
  ki::Vector,
  out::Dict,
  out_field::String,
  out_large=nothing,
  nlevel_large=nothing,
)
  @assert(occursin("t_", out_field))
  tau_r = "$(out["100Float64"]["rlevel"][ri[exp_i]].num)" * string((out["100Float64"]["rlevel"][ri[exp_i]].den == 1 ? "" : "/" * string(out["100Float64"]["rlevel"][ri[exp_i]].den)))
  tau_k = "$(out["100Float64"]["klevel"][ki[exp_i]].num)" * string((out["100Float64"]["klevel"][ki[exp_i]].den == 1 ? "" : "/" * string(out["100Float64"]["klevel"][ki[exp_i]].den)))
  nlevel_b = [@sprintf("%d", floor(log10(nlevelj[j]))) for j in eachindex(nlevelj)]
  data_mn = hcat(
    [
      [
        haskey(out[string((nlevelj[j],DT)...)], algl[l]) ? nanmean(
          out[string((nlevelj[j],DT)...)][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
        ) : +Inf
        for j in eachindex(nlevelj)
      ]
      for l in eachindex(algl)
    ]...
  )'
  data_sd = hcat(
    [
      [
        haskey(out[string((nlevelj[j],DT)...)], algl[l]) ? nanstd(
          out[string((nlevelj[j],DT)...)][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
        ) : +Inf
        for j in eachindex(nlevelj)
      ]
      for l in eachindex(algl)
    ]...
  )'
  if !isnothing(out_large) && !isnothing(nlevel_large)
    data_mn_large = hcat(
      [
        [
          haskey(out_large[string((nlevel_large[j],DT)...)], algl[l]) ? nanmean(
            out_large[string((nlevel_large[j],DT)...)][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
          ) : +Inf
          for j in eachindex(nlevel_large)
        ]
        for l in eachindex(algl)
      ]...
    )'
    data_sd_large = hcat(
      [
        [
          haskey(out_large[string((nlevel_large[j],DT)...)], algl[l]) ? nanstd(
            out_large[string((nlevel_large[j],DT)...)][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
          ) : +Inf
          for j in eachindex(nlevel_large)
        ]
        for l in eachindex(algl)
      ]...
    )'
    large_mask = nlevelj .∈ Ref(nlevel_large)
    data_mn[3:5,large_mask] .= data_mn_large[3:5,:]
    data_sd[3:5,large_mask] .= data_sd_large[3:5,:]
    # println(data_sd[3:5,large_mask][data_mn[3:5,large_mask] .> 10_000])
    @views data_sd[3:5,large_mask][data_mn[3:5,large_mask] .> 10_000] .= +Inf
    println(data_sd)
  end
  algbest = [findall(d .== minimum(d)) for d in eachcol(data_mn)]
  if !isnothing(out_large) && !isnothing(nlevel_large)
    return [
      ["", MultiColumn(3, :l, "Experiment $(exp_i): " * L"\tau_r=%$(tau_r),\,\tau_k^\complement=%$(tau_k)")],
      CMidRule(2, 6),
      ["", [L"n=10^{%$(nlevel_b[j])}" for j in eachindex(nlevelj)]...],
      Rule(:mid),
      [
        [algl_names[l], [(algl[l] in algbest[j] ? L"\textbf{%$(sci_1a(data_mn[l,j]))}" : "$(sci_1a(data_mn[l,j]))") * (nlevelj[j] ∈ nlevel_large && l >= 3 ? "*" : "") * " (" * sci_1b(data_sd[l,j]) * (nlevelj[j] ∈ nlevel_large && l >= 3 ? "*" : "") * ")" for j in eachindex(nlevelj)]...]
        for l in eachindex(algl)
      ]...,
    ]
  else
    return [
      ["", MultiColumn(3, :l, "Experiment $(exp_i): " * L"\tau_r=%$(tau_r),\,\tau_k^\complement=%$(tau_k)")],
      CMidRule(2, 6),
      ["", [L"n=10^{%$(nlevel_b[j])}" for j in eachindex(nlevelj)]...],
      Rule(:mid),
      [
        [algl_names[l], [(algl[l] in algbest[j] ? L"\textbf{%$(sci_1a(data_mn[l,j]))}" : "$(sci_1a(data_mn[l,j]))") * " (" * sci_1b(data_sd[l,j]) * ")" for j in eachindex(nlevelj)]...]
        for l in eachindex(algl)
      ]...,
      # [
      #   ["Algorithm $(algl[l])", [mnstd(out[string((nlevelj[j],DT)...)][algl[l]]["t_total"][:,ri[exp_i],ki[exp_i]]) for j in eachindex(nlevelj)]...]
      #   for l in eachindex(algl)
      # ]...,
      # Rule(:bottom),
    ]
  end
end

function experiment_pct(exp_i::Integer,
  nlevelj::Vector,
  algl::Vector,
  algl_names::Vector,
  ri::Vector,
  ki::Vector,
  out::Dict,
  out_field::String
)
  tau_r = "$(out["100Float64"]["rlevel"][ri[exp_i]].num)" * string((out["100Float64"]["rlevel"][ri[exp_i]].den == 1 ? "" : "/" * string(out["100Float64"]["rlevel"][ri[exp_i]].den)))
  tau_k = "$(out["100Float64"]["klevel"][ki[exp_i]].num)" * string((out["100Float64"]["klevel"][ki[exp_i]].den == 1 ? "" : "/" * string(out["100Float64"]["klevel"][ki[exp_i]].den)))
  nlevel_a = [@sprintf("%d", log2(nlevelj[j])) for j in eachindex(nlevelj)]
  nlevel_b = [@sprintf("%d", floor(log10(nlevelj[j]))) for j in eachindex(nlevelj)]
  data_mn = hcat(
    [
      [
        haskey(out[string((nlevelj[j],DT)...)], algl[l]) ? nanmean(
          out[string((nlevelj[j],DT)...)][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
        ) : +0.0
        for j in eachindex(nlevelj)
      ]
      for l in eachindex(algl)
    ]...
  )'
  j3 = findall(nlevelj .> maxn_grid)
  j4 = findall(nlevelj .> maxn_gurobi)
  j5 = findall(nlevelj .> maxn_gurobi)
  # algbest = [argmax(d) for d in eachcol(data_mn)]
  algbest = [findall(d .== maximum(d)) for d in eachcol(data_mn)]
  data_mn[findfirst(algl.==3),j3] .= NaN
  data_mn[findfirst(algl.==4),j4] .= NaN
  data_mn[findfirst(algl.==5),j5] .= NaN
  return [
    ["", MultiColumn(5, :l, "Experiment $(exp_i): " * L"\tau_r=%$(tau_r),\,\tau_k=%$(tau_k)")],
    CMidRule(2, 6),
    ["", [L"n=10^{%$(nlevel_b[j])}" for j in eachindex(nlevelj)]...],
    Rule(:mid),
    [
      [algl_names[l], [(algl[l] in algbest[j] ? L"\textbf{%$(pct_safe(data_mn[l,j]))}" : "$(pct_safe(data_mn[l,j]))") for j in eachindex(nlevelj)]...]
      for l in eachindex(algl)
    ]...,
    # [
    #   ["Algorithm $(algl[l])", [mnstd(out[string((nlevelj[j],DT)...)][algl[l]]["t_total"][:,ri[exp_i],ki[exp_i]]) for j in eachindex(nlevelj)]...]
    #   for l in eachindex(algl)
    # ]...,
    # Rule(:bottom),
  ]
end

#
# table: time
#

latex_tabular(PROJPATH*"figures/table_time.tex",
  Tabular("llllll"),
  [
    Rule(:top),
    experiment_time(1, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    experiment_time(2, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    experiment_time(3, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    experiment_time(4, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    experiment_time(5, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    experiment_time(6, nlevelj, algl, algl_names, ri, ki, out, "t_total", out_large, nlevel_large)...,
    Rule(:mid),
    sort_time_row(out, nlevelj),
    partial_sort_time_row(out, nlevelj),
    Rule(:bottom),
  ]
)

#
# table: solution quality
#

latex_tabular(PROJPATH*"figures/table_bestfeas.tex",
  Tabular("llllll"),
  [
    Rule(:top),
    experiment_pct(1, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:mid),
    experiment_pct(2, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:mid),
    experiment_pct(3, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:mid),
    experiment_pct(4, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:mid),
    experiment_pct(5, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:mid),
    experiment_pct(6, nlevelj, algl, algl_names, ri, ki, out, "bestfeas")...,
    Rule(:bottom),
  ]
)

latex_tabular(PROJPATH*"figures/table_bestfeasunsorted.tex",
  Tabular("llllll"),
  [
    Rule(:top),
    experiment_pct(1, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:mid),
    experiment_pct(2, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:mid),
    experiment_pct(3, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:mid),
    experiment_pct(4, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:mid),
    experiment_pct(5, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:mid),
    experiment_pct(6, nlevelj, algl, algl_names, ri, ki, out, "bestfeasunsorted")...,
    Rule(:bottom),
  ]
)

#
# figure: relative time
#

fig, ax = subplots(figsize=(6, 6))
rc("font", size=16)

n = 10_000_000
# n = 1_000_000
# n = 100_000
key = string((n,DT)...)
rlevel = out[key]["rlevel"]
klevel = out[key]["klevel"]
nrep, nr, nk = size(out[key][1]["t_total"])

p_21_t_total = (
  # "(" * L"t_{\mathrm{PLCP}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * L"n=" * @sprintf("%0.1e", n),
  "",
  "21_t_total"
)
o21_t_total = reshape(nanmean(out[key][2]["t_total"] ./ out[string((n,DT)...)][1]["t_total"], dims=1), (nr,nk))
o21_t_total[o21_t_total .< 1] .= -1 ./ o21_t_total[o21_t_total .< 1]
plot_grid(o21_t_total, p_21_t_total...)

p_21_t_primal = (
  # "(" * L"t_{\mathrm{PLCP}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * @sprintf("%0.1e", n),
  "",
  "21_t_primal"
)
o21_t_primal = reshape(nanmean(out[key][2]["t_primal"] ./ out[string((n,DT)...)][1]["t_primal"], dims=1), (nr,nk))
o21_t_primal[o21_t_primal .< 1] .= -1 ./ o21_t_primal[o21_t_primal .< 1]
plot_grid(o21_t_primal, p_21_t_primal...)

p_21_t_pivot = (
  # "(" * L"t_{\mathrm{PLCP}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * @sprintf("%0.1e", n),
  "",
  "21_t_pivot"
)
o21_t_pivot = reshape(nanmean(out[key][2]["t_run"] ./ out[string((n,DT)...)][1]["t_run"], dims=1), (nr,nk))
o21_t_pivot[o21_t_pivot .< 1] .= -1 ./ o21_t_pivot[o21_t_pivot .< 1]
plot_grid(o21_t_pivot, p_21_t_pivot...)

n = 1000
n = 100_000
key = string((n,DT)...)
rlevel = out[key]["rlevel"]
klevel = out[key]["klevel"]
nrep, nr, nk = size(out[key][1]["t_total"])
best_grbs = +Inf
best_grbu = +Inf

p_31_t_total = (
  # "(" * L"t_{\mathrm{GRID}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * L"n=" * @sprintf("%0.1e", n),
  "",
  "31_t_total"
)
o31_t_total = reshape(nanmean(out[key][3]["t_total"] ./ out[string((n,DT)...)][1]["t_total"], dims=1), (nr,nk))
o31_t_total[o31_t_total .< 1] .= -1 ./ o31_t_total[o31_t_total .< 1]
plot_grid(o31_t_total, p_31_t_total...)

p_41_t_total = (
  # "(" * L"t_{\mathrm{GRBS}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * L"n=" * @sprintf("%0.1e", n),
  "",
  "41_t_total"
)
o41_t_total = reshape(nanmean(out[key][4]["t_run"] ./ out[string((n,DT)...)][1]["t_total"], dims=1), (nr,nk))
o41_t_total[o41_t_total .< 1] .= -1 ./ o41_t_total[o41_t_total .< 1]
best_grbs = min(best_grbs, minimum(o41_t_total))
plot_grid(o41_t_total, p_41_t_total...)

p_51_t_total = (
  # "(" * L"t_{\mathrm{GRBU}} / t_{\mathrm{ESGS}}" * ")x speedup @ " * L"n=" * @sprintf("%0.1e", n),
  "",
  "51_t_total"
)
o51_t_total = reshape(nanmean(out[key][5]["t_run"] ./ out[string((n,DT)...)][1]["t_total"], dims=1), (nr,nk))
o51_t_total[o51_t_total .< 1] .= -1 ./ o51_t_total[o51_t_total .< 1]
best_grbu = min(best_grbu, minimum(o51_t_total))
plot_grid(o51_t_total, p_51_t_total...)

close(fig)

#
# figure: time vs n
#

ri = [6,8,11,6,8,11]; out["rlevel"][ri]
ki = [2,2,2,4,4,4]; out["klevel"][ki]
nlevel_all = 10 .^ collect(1:7)
nlarge = [10^6, 10^7]

# fig, ax = subplots(2, 2, figsize=(15, 15))
fig, ax = subplots(figsize=(8, 4)) # wide
# fig, ax = subplots(figsize=(4, 6)) # tall
rc("font", size=18)

for exp_i in eachindex(ri)
  plotname = "exp$(exp_i)_t_vs_n"
  # plot_loglog_r_k(exp_i, plotname, "t_total", ri, ki, algl, algl_names, nlevelj)
  plot_loglog_r_k(exp_i, plotname, "t_total", ri, ki, algl, algl_names, nlevel_all, out_large, nlarge)
end

# save legend
ax.axis("off")
legend = fig.legend(ncol=6, loc="lower right", fontsize=18)
legfig = legend.figure
# legfig.canvas.draw()
# bbox  = legend.get_window_extent().transformed(legfig.dpi_scale_trans.inverted())
# legfig.savefig(PROJPATH * "/figures/$(plotname)_legend.pdf", bbox_inches=bbox, pad_inches=0.1)
legfig.savefig(PROJPATH * "/figures/$(plotname)_legend.pdf", bbox_inches="tight", pad_inches=0.1)

close(fig)

#
# worst t2/t1
#

worst = +Inf
for n in nlevel
  key = string((n,DT)...)
  o21_t_total = reshape(nanmean(out[key][2]["t_total"] ./ out[key][1]["t_total"], dims=1), (nr,nk))
  worst = min(worst, minimum(o21_t_total))
end
1-1/worst

#
# averaged worst t3 / (slow(t1, t2)) for tau in (0, 0.1) and r in (0.1, 0.999)
#

worst_ratio_3 = +Inf
worst_ratio_45 = +Inf
best_ratio_3 = -Inf
best_ratio_45 = -Inf
expers = [string(n)*"Float64" for n in 10 .^ collect(1:5)]
d45_d12 = zeros(12, 5)
d3_d12 = zeros(12, 5)
for exper in expers
  for i in 1:100
    d5 = readdlm(DATAPATH * exper * "/out_5_t_total_$i.csv")[1:end,1:5]
    d4 = readdlm(DATAPATH * exper * "/out_4_t_total_$i.csv")[1:end,1:5]
    d3 = readdlm(DATAPATH * exper * "/out_3_t_total_$i.csv")[1:end,1:5]
    d2 = readdlm(DATAPATH * exper * "/out_2_t_total_$i.csv")[1:end,1:5]
    d1 = readdlm(DATAPATH * exper * "/out_1_t_total_$i.csv")[1:end,1:5]
    d3_d12 .+= d3 ./ max.(d1, d2)
    d45_d12 .+= min.(d4,d5) ./ max.(d1, d2)
  end
  d3_d12 ./= 100
  d45_d12 ./= 100
  worst_ratio_exper_3 = minimum(d3_d12)
  if minimum(d3_d12) < 1
    println(exper)
  end
  worst_ratio_exper_45 = minimum(d45_d12)
  worst_ratio_3 = min(worst_ratio_3, worst_ratio_exper_3)
  worst_ratio_45 = min(worst_ratio_45, worst_ratio_exper_45)
  best_ratio_exper_3 = maximum(d3_d12)
  best_ratio_exper_45 = maximum(d45_d12)
  best_ratio_3 = max(best_ratio_3, best_ratio_exper_3)
  best_ratio_45 = max(best_ratio_45, best_ratio_exper_45)
end
worst_ratio_3
worst_ratio_45
best_ratio_3
best_ratio_45
# julia> worst_ratio_3 #! ok, noisy for n <= 100
# 0.7586453137125331
# julia> worst_ratio_45
# 132.98182996259732
# julia> best_ratio_3
# 57737.58782241328
# julia> best_ratio_45
# 430011.30424307287

worst_ratio_3 = +Inf
worst_ratio_45 = +Inf
best_ratio_3 = -Inf
best_ratio_45 = -Inf
expers = [string(n)*"Float64" for n in 10 .^ collect(5:5)]
d45_d12 = zeros(12, 5)
d3_d12 = zeros(12, 5)
for exper in expers
  for i in 1:100
    d5 = readdlm(DATAPATH * exper * "/out_5_t_total_$i.csv")[1:end,1:5]
    d4 = readdlm(DATAPATH * exper * "/out_4_t_total_$i.csv")[1:end,1:5]
    d3 = readdlm(DATAPATH * exper * "/out_3_t_total_$i.csv")[1:end,1:5]
    d2 = readdlm(DATAPATH * exper * "/out_2_t_total_$i.csv")[1:end,1:5]
    d1 = readdlm(DATAPATH * exper * "/out_1_t_total_$i.csv")[1:end,1:5]
    d3_d12 .+= d3 ./ max.(d1, d2)
    d45_d12 .+= min.(d4,d5) ./ max.(d1, d2)
  end
  d3_d12 ./= 100
  d45_d12 ./= 100
  worst_ratio_exper_3 = minimum(d3_d12)
  if minimum(d3_d12) < 1
    println(exper)
  end
  worst_ratio_exper_45 = minimum(d45_d12)
  worst_ratio_3 = min(worst_ratio_3, worst_ratio_exper_3)
  worst_ratio_45 = min(worst_ratio_45, worst_ratio_exper_45)
  best_ratio_exper_3 = maximum(d3_d12)
  best_ratio_exper_45 = maximum(d45_d12)
  best_ratio_3 = max(best_ratio_3, best_ratio_exper_3)
  best_ratio_45 = max(best_ratio_45, best_ratio_exper_45)
end
worst_ratio_3
worst_ratio_45
best_ratio_3
best_ratio_45

#=
...noisy
#
# instance-wise worst t3 / (slow(t1, t2)) for tau in (0, 0.1) and r in (0.1, 0.999)
#
worst_ratio = +Inf
worst_ratio_exper = +Inf
expers = [string(n)*"Float64" for n in 10 .^ collect(3:3)]
for exper in expers
  for i in 1:100
    d3 = readdlm(DATAPATH * exper * "/out_3_t_total_$i.csv")[1:end,1:5]
    d2 = readdlm(DATAPATH * exper * "/out_2_t_total_$i.csv")[1:end,1:5]
    d1 = readdlm(DATAPATH * exper * "/out_1_t_total_$i.csv")[1:end,1:5]
    worst_ratio_exper = minimum(d3./max.(d2,d1))
    println(worst_ratio_exper)
    worst_ratio = min(worst_ratio, worst_ratio_exper)
  end
  println(worst_ratio)
  println(exper)
  println()
end
worst_ratio
=#