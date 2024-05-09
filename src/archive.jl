
ODEparams = [1000.0, 2.0, 200.0, 2.5]
params = vcat(ones(nres+2), ODEparams)
times = 1.0:Δt:1095.0
Q = gr4j(params, times)


# loss_function(p) = NSE_loss(gr4j, p, train_Y, train_points)
# loss_function(params)

# initial_params = params


opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

optimizer = ADAM(0.1)
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=50)

# #test data
# out_params = sol.u
# times = train_ + 1:Δt:data_points[end]

# Q_nn = gr4j(out_params, times)

# plot(times, Q_nn)
# plot!(times, df[(train_ +1):end,:streamflow], dpi = 300)
# plot!(times,Q_nn)



# savefig("chepe_plot_gr4j_itr200.png")


ODEparams = [1000.0, 2.0, 200.0, 2.5, 10.0, 13.0, 15.0, 0.0]
params = vcat(ones(nres+3), ODEparams)
times = 1.0:Δt:data_points[end]
Qsnow, sol = gr4jsnow(params, times)


# loss_function(p) = NSE_loss(gr4jsnow, p, train_Y, train_points)
# loss_function(params)

# initial_params = params


# opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

# opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

# optimizer = ADAM(0.1)
# sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=70)

# #test data
# out_params = sol.u
# times = train_ + 1:Δt:data_points[end]

# Q_nn = gr4jsnow(out_params, times)

# plot(times, Q_nn)
# plot!(times, df[(train_ +1):end,:streamflow], dpi = 300)
# plot!(times,Q_nn)



# savefig("chepe_plot_gr4jsnow_itr200.png")






ann, initial_NN_params = initializeNN()
ODEparams = [1000.0, 2.0, 200.0, 2.5]
ODE_states = ones(nres+3)
# initial_params = (ODE_states, ODEparams, initial_NN_params)
initial_params = ComponentArray(ODE_states=ODE_states, ODEparams=ODEparams, NNparams=initial_NN_params)


#Load initial_params from JLD:
initial_params = load_object("chepe_params_tem.jld")

gr4jsnowNN(initial_params, train_points, ann)
UDE_model(params, output_times) = gr4jsnowNN(params, output_times, ann)

loss_function(p) = NSE_loss(UDE_model, p, train_Y, train_points)
loss_function(initial_params)

opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

# optimizer = LBFGS()
optimizer = ADAM(0.001)
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=100000)

# out_params = sol.u
# times = 1.0:Δt:train_

# Q_nn = gr4jsnowNN(out_params, times, ann)
# plot(times, Q_nn)
# plot!(times, df[1:train_,:streamflow], dpi = 300)
# plot!(times,Q_nn)

#test data
times = train_ + 1:Δt:data_points[end]

Q_nn = gr4jsnowNN(out_params, times, ann)

plot(times, Q_nn)
plot!(times, df[(train_ +1):end,:streamflow], dpi = 300)
plot!(times,Q_nn)




# savefig("chepe_plot_gr4jsnowNN_itr70.png")
