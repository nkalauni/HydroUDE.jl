include("startup.jl")
include("load_data.jl")
include("utils.jl")
include("models.jl")
include(pwd()*"/src/optimize.jl")

#=========================================================================#
# Vanilla GR4J
# Initialize 4 params + initial states.
ODEparams = [1000.0, 2.0, 200.0, 2.5]
ODEstates = ones(nres+2)

initial_params = vcat(ODEstates, ODEparams)
Wrapper_model(p,t) = GR4J_model(p,t)
Wrapper_model(initial_params, train_points)
# optm_parameters = optimize_model(Wrapper_model, initial_params, maxitr=2)
optm_parameters = load_object("optim_vars/gr4j_state-params.jld")

#save and plot the model
save_model(optm_parameters, "t")
plot_model(Wrapper_model, optm_parameters, portion="train", name="gr4j_sp")
plot_model(Wrapper_model, optm_parameters, portion="test", name="gr4j_sp")

#=========================================================================#
# * optimize 4 params starting with optimized initial states.
ODEstates = optm_parameters[1:nres+2]    # Read ODEstates from previous optimization process

Wrapper_model(p, t) = GR4J_model(p, t, ODEstates)
# gr4j_params = optimize_model(Wrapper_model, ODEparams, maxitr = 2)
gr4j_params = load_object("optim_vars/gr4j_params.jld")

#save and plot the model
save_model(gr4j_params, "t")
plot_model(Wrapper_model, gr4j_params, portion="train", name="gr4j_p")
plot_model(Wrapper_model, gr4j_params, portion="test", name="gr4j_p")

#=========================================================================#
#=========================================================================#
# Vanilla GR4Jsnow
# Initialize 8 params + initial states.
ODEparams = [1000.0, 2.0, 200.0, 2.5, 10.0, 13.0, 15.0, 0.0]
ODEstates = ones(nres+3)

initial_params = vcat(ODEstates, ODEparams)
Wrapper_model(p,t) = gr4jsnow(p,t)
Wrapper_model(initial_params, train_points)
# gr4jsnow_sparams = optimize_model(Wrapper_model, initial_params, maxitr=2)
gr4jsnow_sparams = load_object("optim_vars/gr4jsnow_state-params.jld")

#save and plot the model
save_model(gr4jsnow_sparams, "t")
plot_model(Wrapper_model, gr4jsnow_sparams, portion="train", name="gr4jsnow_sp")
plot_model(Wrapper_model, gr4jsnow_sparams, portion="test", name="gr4jsnow_sp")

#=========================================================================#
# * optimize 8 params starting with optimized initial states from GR4Jsnow.

ODEstates = gr4jsnow_sparams[1:nres+3]   # Read ODEstates from previous optimization process

Wrapper_model(p, t) = gr4jsnow(p, t, ODEstates)
Wrapper_model(ODEparams, train_points)
# gr4jsnow_params = optimize_model(Wrapper_model, ODEparams, maxitr = 2)
gr4jsnow_params = load_object("optim_vars/gr4jsnow_params.jld")
#save and plot the model
save_model(gr4jsnow_params, "t")

Wrapper_model(gr4jsnow_params, train_points)

portion = "train"
p = gr4jsnow_params
model = Wrapper_model
Q_nn = model(p, test_points)
plot(test_points, test_Y, dpi = 300)
plot!(test_points, Q_nn, title=portion*"_NSE: "*string(-NSE_loss(model, p, test_Y, test_points)))
plot!(test_points, Q_nn)
savefig("plots/t.png")

plot_model(Wrapper_model, gr4jsnow_params, portion="train", name="gr4jsnow_p")
plot_model(Wrapper_model, gr4jsnow_params, portion="test", name="gr4jsnow_p")

#=========================================================================#
#=========================================================================#
# ## GR4JsnowNN--
# * optimize NN params starting with optimized initial states from GR4Jsnow.

ann, initial_NN_params = initializeNN()
ODEparams = [1000.0, 2.0, 200.0, 2.5] # redefine ODE params
ODEstates = gr4jsnow_sparams[1:nres+3]  #use ODE states from gr4jsnow model.

# ODEstates = ones(nres+3)
# initial_sparams = ComponentArray(ODEstates=ODEstates, ODEparams=ODEparams, NNparams=initial_NN_params)
initial_params = ComponentArray(ODEparams = ODEparams, NNparams = initial_NN_params)

Wrapper_model(p, t) = gr4jsnowNN(p, t, ann, ODEstates)
Wrapper_model(initial_params, train_points)
optm_parameters = optimize_model(Wrapper_model, initial_params, maxitr = 2)
# optm_parameters = load_object("optim_vars/gr4jsnowNN_params.jld")

#save and plot the model
save_model(optm_parameters, "t")
plot_model(Wrapper_model, optm_parameters, portion="train", name="gr4jsnowNN_p")
plot_model(Wrapper_model, optm_parameters, portion="test", name="gr4jsnowNN_p")
