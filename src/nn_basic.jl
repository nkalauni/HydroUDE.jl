using CSV
using DataFrames, Dates, Statistics
using Revise

using Interpolations


df = CSV.read("data/chepe_data.csv", DataFrame; dateformat="mm/dd/yyyy")
data_points = collect(1:nrow(df))
 
interpolate_method = SteffenMonotonicInterpolation()
 
temp = interpolate(data_points, df[:,:tavg], interpolate_method)
prec = interpolate(data_points, df[:,:precipitation], interpolate_method)
pet = interpolate(data_points, df[:,:pet], interpolate_method)

println(pwd())
println("0")
