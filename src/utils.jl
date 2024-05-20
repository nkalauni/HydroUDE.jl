
# generic loss function:
function NSE_loss(model, params, target_data, target_time)

    predicted_data = model(params, target_time)
    NSE(obs, pred) = 1 - sum((obs .- pred).^2) / sum((obs .- mean(obs)).^2)
    if size(predicted_data) == size(target_data)
        # println("size milyo!!! Yesss..... Data_size = $(size(target_data))......ODE solution size = $(size(predicted_data))")
        loss = -NSE(target_data, predicted_data)
    else
        println("Size milena re ODE solve garda....Data_size = $(size(target_data))......ODE solution size = $(size(predicted_data))")
        loss = Inf
    end

    return loss

end

function callback(model, p, l, losses)

    println("Updating plot")
    #txt save
    open("runtime_env/runtime_params.txt", "w") do file
        write(file, "NSE: "*string(-l)*", "*string(p.u))
    end
    
    #jld save
    save_object("runtime_env/runtime_params.jld", p.u)

    #iteration counter
    push!(losses, -l)

    #image save
    Q_nn = model(p.u, test_points)
    plot(test_points, test_Y, dpi = 300)
    plot!(test_points, Q_nn, title="NSE[Itr: $(length(losses))]: "*string(-l))
    savefig("runtime_env/runtime_plot.png")

    #Pring log in terminal
    println("NSE[Itr: $(length(losses))]: "*string(-l))
    return false

end
