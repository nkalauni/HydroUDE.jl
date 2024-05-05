includet("startup.jl")
include("load_data.jl")
includet("utils.jl")
includet("models.jl")

# ## GR4J
# optimize 4 params + initial states.
ODEparams = [1000.0, 2.0, 200.0, 2.5]
ODEstates = ones(nres+2)
initial_params = vcat(ODEstates, ODEparams)

# define loss function
loss_function(p) = NSE_loss(gr4j, p, train_Y, train_points)

# define optimization  function
opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

optimizer = ADAM(0.1)
# sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=50)
#save optimized parameters
# save_object("optim_vars/gr4j_state-params.jld", sol.u)



# * optimize 4 params starting with optimized initial states.
parameters = load_object("optim_vars/gr4j_state-params.jld")
ODEstates = parameters[1:nres+2]
ODEparams = parameters[nres+3:end]

GR4J_model(ODEparams, train_points, ODEstates)


#Model wrapper
Wrapper_model(p, t) = GR4J_model(p, t, ODEstates)

Wrapper_model(ODEparams, train_points)
# define loss function
loss_function(p) = NSE_loss(Wrapper_model, p, train_Y, train_points)
loss_function(ODEparams)

callback_function(p,l) = callback(Wrapper_model, p, l)

# define optimization  function
opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

opt_problem = Optimization.OptimizationProblem(opt_func, ODEparams)

optimizer = ADAM(0.1)
sol = Optimization.solve(opt_problem, optimizer, callback = callback_function, maxiters=1000)

#save optimized parameters
# save_object("optim_vars/gr4j_params.jld", load_object("chepe_params.jld"))

optm_params = load_object("optim_vars/gr4j_params.jld")

div = Int(round(length(train_points)/2))
train_points = train_points[1:div]
train_Y = train_Y[1:div]

Q_nn = Wrapper_model(optm_params, train_points)
plot(train_points, train_Y, dpi = 300)
plot!(train_points, Q_nn, title="NSE: "*string(-NSE_loss(Wrapper_model, optm_params, train_Y, train_points)))



Q_nn = Wrapper_model(optm_params, test_times)
plot(test_times, df[(train_ +1):end,:streamflow], dpi = 300)
plot!(test_times,Q_nn, title="NSE: "*string(-NSE_loss(Wrapper_model, optm_params, test_Y, test_points)))







# ## GR4JsnowNN
# * optimize NN params starting with optimized 4 params from GR4J.
