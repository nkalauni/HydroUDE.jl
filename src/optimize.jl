function optimize_model(Wrapper_model, initial_params; optimizer = ADAM(0.1), maxitr = 20, save = false)

    # define loss function
    loss_function(p) = NSE_loss(Wrapper_model, p, train_Y, train_points)

    # define callack function
    callback_function(p,l) = callback(Wrapper_model, p, l)

    # define optimization  function
    opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

    # formulate optimization problem
    opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

    sol = Optimization.solve(opt_problem, optimizer, callback=callback_function, maxiters=maxitr)
    return sol.u
    
    
    # save optimized parameters
    # save_object("optim_vars/gr4j_state-params.jld", sol.u)
    # save_object("optim_vars/gr4j_state-params.jld", load_object("chepe_params.jld"))
    # parameters = load_object("optim_vars/gr4j_state-params.jld")

end

function save_model(model_params, filename)
    path = pwd()*"/optim_vars/"*filename*".jld"
    save_object(path, model_params)
end

function plot_model(model, p; portion="test", save= true, name="noname")
    if portion == "test"
        Q_nn = model(p, test_points)
        plot(test_points, test_Y, dpi = 300)
        plot!(test_points, Q_nn, title="NSE: "*string(-NSE_loss(model, p, test_Y, test_points)))
    else
        Q_nn = model(p, train_points)
        plot(train_points, train_Y, dpi = 300)
        plot!(train_points, Q_nn, title="NSE: "*string(-NSE_loss(model, p, train_Y, train_points)))
    end

    if save == true
        savefig("plots/"*name*"_model_"*portion*"_plot.png")
    end
    
end

function plot_model(model, p, points, values)
    Q_nn = model(optm_params, points)
    plot(points, value, dpi = 300)
    plot!(points, Q_nn, title="NSE: "*string(-NSE_loss(model, p, values, points)))
end
