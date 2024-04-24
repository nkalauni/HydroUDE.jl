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

df

NSE(pred, obs) = 1 - sum((pred .- obs).^2) / sum((obs .- mean(obs)).^2)

function NSE_loss(pred_model, params, batch, time_batch)

    pred, = pred_model(params, time_batch)
    loss = -NSE(pred,batch)

    return loss, pred
end


input_var = ["precipitation","tmin","tmax","tavg","wind","pet"]
output_var = "streamflow"
