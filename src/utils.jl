# generic loss function:
function NSE_loss(model, params, target_data, target_time)

    predicted_data = model(params, target_time)
    NSE(obs, pred) = 1 - sum((obs .- pred).^2) / sum((obs .- mean(obs)).^2)
    loss = -NSE(target_data, predicted_data)

    return loss 

end

function callback(model, p, l)

    #txt save
    open("chepe_log_file.txt", "w") do file
        write(file, "NSE: "*string(-l)*", "*string(p.u))
    end

    #json save
    # open("chepe_.json", "w") do f
    #     JSON.print(f, p.u, 4)  # The '4' here specifies the indentation level for pretty printing
    # end

    #jld save
    save_object("chepe_params.jld", p.u)

    #image save
    println("Updating plot")
    times = train_ + 1:Δt:data_points[end]
    
    # Q_nn = gr4jsnowNN(p.u, times, ann)
    Q_nn = model(p.u, times)

    
    plot(times, df[(train_ +1):end,:streamflow], dpi = 300)
    plot!(times,Q_nn, title="NSE: "*string(-l))
    
    savefig("src/runtime_plot.png")

    println("NSE: "*string(-l))
    return false

end
