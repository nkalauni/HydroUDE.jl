
# Read data and preprocessing
df = CSV.read("data/modified_chepe_data.csv", DataFrame; dateformat="mm/dd/yyyy")

# Interpolate to generate continuous data
data_points = collect(1:nrow(df))
interpolate_method = SteffenMonotonicInterpolation()

temp = interpolate(data_points, df[:,:temperature], interpolate_method)
prec = interpolate(data_points, df[:,:precipitation], interpolate_method)
pet = interpolate(data_points, df[:,:pet], interpolate_method)
flow = df[:, :streamflow]

# Split data for training and testing 
const split_ratio = 0.7
train_ = Int(round(split_ratio * nrow(df)))
test_ = nrow(df) - train_

# Split training data into feature and target
train_Y = flow[1:train_]
train_points = collect(1:train_)

test_Y = flow[(train_ +1):end]
test_points = collect((train_ +1):nrow(df))


# Define time interval for ODE solver and redefine data_points accordingly
const Δt = 1.0
data_points = collect(1:Δt:nrow(df))
train_times = 1:Δt:train_
test_times = train_ + 1:Δt:data_points[end]
