using CSV
using DataFrames
using Plots
using DelimitedFiles

dat = CSV.read("output.csv", DataFrame)
p = plot(dat.time, dat.data, seriestype = :line, title = "Corr Function", xlabel = "log(t)", ylabel = "log|C(t)|")
savefig(p, "test.png")

# data = readdlm("time_correlation_data.txt", ' ')  # Replace ' ' with '\t' for tab-delimited, or ',' for comma-delimited
# x = data[1:end, 1] 
# y = data[1:end, 2]
# y = abs.(y)
# plot(x, y, seriestype = :scatter, title = "Scatter Plot", xlabel = "X Axis", ylabel = "Y Axis")
# savefig("plot_from_txt.png")
