using PyCall, PyPlot
using NaNStatistics
# plt.switch_backend("agg") #!no
@pyimport mpl_toolkits.axes_grid1 as ag1
rc("text", usetex=true)
rc("font", family="serif")
rc("mathtext", fontset="cm")
rc("legend", fancybox=false)
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

sci_1c(d) = begin
  """one-exponent-digit and zero-base-decimaldigit if d>10; o.w. if 1<d<10 use d"""
  expo = floor(log10(abs(d)))
  base = d * 10^-expo
  if d<=10
    @sprintf("%0.1f", d)
  elseif isnan(base) || isinf(base) || isnan(expo) || isinf(expo)
    L"-"
  else
    @sprintf("%0.1fe%d", base, expo)
  end
end
function plot_grid(data::Matrix, plottitle::String, plotname::String)
  ax = subplot(1,1,1)
  ax.cla()
  divider = ag1.axes_divider.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  # cmap = get_cmap("viridis")
  cmap = get_cmap("Reds")
  cmap.set_under("black")
  if occursin("/", plottitle)
    # ratio
    im = ax.imshow(data, cmap=cmap, vmin=1.0) # https://stackoverflow.com/a/63482254/6272122
  else
    # non-ratio
    # im = ax.imshow(data, cmap=cmap)
    im = ax.imshow(data, cmap=cmap, vmin=1.0) # https://stackoverflow.com/a/63482254/6272122
  end
  ax.set_ylabel(L"\tau_r,\; (r = \tau_r\cdot \mathsf{M}_{(k)}(x^0))", fontsize=16)
  ax.set_yticks(collect(0:nr-1),rlevel)
  ax.set_xlabel(L"\tau_k^{\mathsf{c}},\; (k = \tau_k^{\mathsf{c}}\cdot n)", fontsize=16)
  ax.set_xticks(collect(0:nk-1), klevel, rotation=45)
  ax.set_title(plottitle)
  cb = fig.colorbar(im, cax=cax, orientation="vertical", extend="min", extendrect=true)
  cb.ax.set_yticklabels([sci_1c(x) for x in cb.get_ticks()])
  fig.tight_layout()
  fig.savefig(PROJPATH * "/figures/$(key...)_$(plotname).pdf")
  if isa(cb,PyObject)
    cb.remove()
  end
end
#=
function plot_grid(data::Matrix, plottitle::String, plotname::String)
  ax = subplot(1,1,1)
  ax.cla()
  divider = ag1.axes_divider.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  # cmap = get_cmap("viridis")
  # cmap.set_under("black",alpha=0.75)
  
  cmap = get_cmap("coolwarm")
  if minimum(data) < 0
    divnorm = matplotlib.colors.TwoSlopeNorm(vcenter=0.0, vmax=maximum(data)+0.01, vmin=minimum(data)-0.01) # https://stackoverflow.com/a/56699813/6272122
  else
    # divnorm = matplotlib.colors.TwoSlopeNorm(vmax=maximum(data)+0.01, vmin=0.0) # https://stackoverflow.com/a/56699813/6272122

    py"""
    import numpy as np
    import matplotlib as mpl
    class MidpointNormalizeFair(mpl.colors.Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
            # I'm ignoring masked values and all kinds of edge cases to make a
            # simple example...

            result, is_scalar = self.process_value(value)
            self.autoscale_None(result)

            vlargest = max( abs( self.vmax - self.midpoint ), abs( self.vmin - self.midpoint ) )
            x, y = [ self.midpoint - vlargest, self.midpoint, self.midpoint + vlargest], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y))
    """
    divnorm = py"""MidpointNormalizeFair( midpoint = 0 )"""
  end
  
  # if occursin("/", plottitle)
    # ratio
    # im = ax.imshow(data, cmap=cmap, vmin=1.0) # https://stackoverflow.com/a/63482254/6272122
    im = ax.imshow(data, cmap=cmap, norm=divnorm) # https://stackoverflow.com/a/63482254/6272122
  # else
    # non-ratio
    # im = ax.imshow(data, cmap=cmap)
  # end
  ax.set_ylabel(L"\tau_r,\; (r = \tau_r\cdot \mathsf{M}_{(k)}(x^0))", fontsize=16)
  ax.set_yticks(collect(0:nr-1),rlevel)
  ax.set_xlabel(L"\tau_k,\; (k = \tau_k\cdot n)", fontsize=16)
  # ax.set_xticks(collect(0:nk-1), klevel, rotation=45)
  ax.set_xticks(collect(0:nk-1), klevel)
  ax.tick_params(axis="x", labelrotation = 45)
  ax.set_title(plottitle)
  cb = fig.colorbar(im, cax=cax, orientation="vertical", extend="min", extendrect=true)
  # cb.ax.set_yticklabels([@sprintf("%0.1e", x) for x in cb.get_ticks()])
  # cb.ax.locator_params(nbins=10)
  # tlabs = round.(collect(range(minimum(data), stop=maximum(data), length=10)), digits=2)
  # tlocs = round.([x[1] for x in collect(range(divnorm(minimum(data)), stop=divnorm(maximum(data)), length=10))], digits=2)
  tlocs = round.([x[1] for x in collect(range(0, stop=1, length=10))], digits=2)
  tlocs = divnorm.inverse(tlocs)
  if minimum(data) < 0
    tlocs[tlocs.<0] .= round.(tlocs[tlocs.<0], digits=1, RoundUp)
  end
  tlocs[tlocs.>0] .= round.(tlocs[tlocs.>0], digits=1, RoundDown)
  cb.set_ticks(tlocs, labels=[sci_1c(x) for x in round.(tlocs, digits=1)])
  # cb.set_ticks(tlocs, [sci_1c(t) for t in tlocs])
  # cb.set_ticks(log.(ts), [sci_1c(t) for t in ts])
  # cb.ax.set_yticklabels([sci_1c(x) for x in cb.get_ticks()])
  # lo = floor(round(minimum(data[.!isnan.(data)]), digits=1))
  # hi = ceil(round(maximum(data[.!isnan.(data)]), digits=-1))
  # cb.set_ticks(sort([collect(range(lo, stop=hi, length=5)); 1]))
  # cb.set_ticklabels(sort([round.(collect(range(lo, stop=hi, length=5)), digits=0); 1.0]))
  fig.tight_layout()
  fig.savefig(PROJPATH * "/figures/$(key...)_$(plotname).pdf")
  if isa(cb,PyObject)
    cb.remove()
  end
end
=#
function plot_loglog_r_k(
  exp_i::Integer,
  plotname::String,
  out_field::String,
  ri::Vector, ki::Vector, algl::Vector, algl_names::Vector, nlevel::Vector,
  out_large=nothing,
  nlarge=nothing,
)
  keyj = [string(n)*string(DT) for n in nlevel]
  ax = subplot(1,1,1)
  ax.cla()
  
  # if fix_type == "r"
  #   data = [
  #     [
  #       [nanmean(out[keyj[j]][algl[l]][out_field][:,ri[exp_i],ki[i]]) for j in eachindex(nlevel)]
  #       for l in eachindex(algl)
  #     ]
  #     for i in eachindex(ki)
  #   ]
  #   for i in eachindex(ki)
  #     for l in eachindex(algl)
  #       data[i][l][isinf.(data[i][l])] .= NaN
  #     end
  #   end
  # else
  #   data = [
  #     [
  #       [nanmean(out[keyj[j]][algl[l]][out_field][:,ri[i],ki[exp_i]]) for j in eachindex(nlevel)]
  #       for l in eachindex(algl)
  #     ]
  #     for i in eachindex(ri)
  #   ]
  #   for i in eachindex(ri)
  #     for l in eachindex(algl)
  #       data[i][l][isinf.(data[i][l])] .= NaN
  #     end
  #   end
  # end
  data = [
    [
      haskey(out[keyj[j]], algl[l]) ? nanmean(
        out[keyj[j]][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
      ) : NaN
      for j in eachindex(nlevel)
    ]
    for l in eachindex(algl)
  ]
  if !isnothing(out_large) && !isnothing(nlarge)
    large_mask = nlevel .∈ Ref(nlarge)
    keyj_large = keyj[large_mask]
    data_large = [
      [
        haskey(out_large[keyj_large[j]], algl[l]) ? nanmean(
          out_large[keyj_large[j]][algl[l]][out_field][:,ri[exp_i],ki[exp_i]]
        ) : NaN
        for j in eachindex(nlarge)
      ]
      for l in eachindex(algl)
    ]
    for alg in 3:5
      data[alg][large_mask] .= data_large[alg]
    end
  end
  for l in eachindex(algl)
    data[l][isinf.(data[l])] .= NaN
  end
  idx = out["sort_time"][:,1] .∈ Ref(nlevel)
  sort_t = vec(nanmean(out["sort_time"][idx, 2:end], dims=2))

  # for i in eachindex(data)
  lines = []
  for l in eachindex(algl)
    linel = ax.loglog(nlevel, data[l], label=algl_names[l], color=cycle[l], marker="o")
    push!(lines, linel)
    # if i == 1
    #   ax.loglog(nlevel, data[i][l], label=algl_names[l], color=cycle[l], marker="o")
    # else
    #   ax.loglog(nlevel, data[i][l], color=cycle[l], marker="o")
    # end
    # lastidx = findlast(.!isnan.(data[i][l]))
    # if fix_type == "r"
    #   v = out["klevel"][ki[i]]
    #   ax.annotate(xy=(nlevel[lastidx], data[i][l][lastidx]), xytext=(1,0), textcoords="offset points", text=L"\tau_k=%$(v)", va="center")
    # else
    #   v = out["rlevel"][ri[i]]
    #   ax.annotate(xy=(nlevel[lastidx], data[i][l][lastidx]), xytext=(1,0), textcoords="offset points", text=L"\tau_r=%$(v)", va="center")
    # end
  end
  # end
  
  ax.loglog(nlevel, sort_t, label="Sort", color="black", linestyle=":")
  ax.grid(true)
  ax.set_xlabel(L"n", fontsize=18)
  ax.set_ylabel("time (sec)", fontsize=18)
  vk = replace(string(out["klevel"][ki[exp_i]]),"//"=>"/")
  vr = replace(string(out["rlevel"][ri[exp_i]]),"//"=>"/")
  ax.set_title("Experiment $exp_i: " * L"\tau_r=%$(vr),\,\tau_k^{\mathsf{c}}=%$(vk)")
  fig.tight_layout()
  fig.savefig(PROJPATH * "/figures/$(plotname).pdf", bbox_inches = "tight", pad_inches=0.02)
end
