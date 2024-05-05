# Load all external packages
cd(@__DIR__)
cd("..")

using Pkg; Pkg.activate(".")

using Revise

using DataFrames, Dates, Statistics
using CSV, DelimitedFiles

# using DifferentialEquations
using OrdinaryDiffEq
using Lux, DiffEqFlux
using LuxCUDA
# using Lux
using ComponentArrays
using SciMLSensitivity

using Optimization
using OptimizationBBO
using OptimizationOptimJL
using OptimizationOptimisers
using OptimizationCMAEvolutionStrategy

using ForwardDiff
using Zygote

using Interpolations

using Plots

using Random
Random.seed!(300)

include("src/helpers/data_loader.jl")

basin_id = "01013500"
data_path = "data/reduced_camels"

df = load_camels_data(lpad(string(basin_id), 8, "0"), data_path)

# df = load_nepal_data("data/chepe_data.csv")

interpolate_method = SteffenMonotonicInterpolation()

const Δt = 1.0
data_points = collect(1:nrow(df))

temp = interpolate(data_points, df[:,:temperature], interpolate_method)
# temp = interpolate(data_points, df[:,"Tmean(C)"], interpolate_method)
prec = interpolate(data_points, df[:,:precipitation], interpolate_method)
# prec = interpolate(data_points, df[:,"Prec(mm/day)"], interpolate_method)
pet = interpolate(data_points, df[:,:pet], interpolate_method)
# pet_arr = PET.(df[:,"Tmean(C)"], df[:, "Daylight(h)"])
# pet = interpolate(data_points, pet_arr, interpolate_method)
flow = df[:, :streamflow]
# flow = df[:, "Flow(mm/s)"]

split_ratio = 0.5
train_ = Int(round(split_ratio * nrow(df)))
test_ = nrow(df) - train_

train_Y = flow[1:train_]
train_points = collect(1:train_)

const α = 2
const β = 5
const γ = 5
const ω = 3.5
const ϵ = 1.5
const Φ = 0.9
const ν = 4/9
const Uₜ = 1
const nres = 11
const SMOOTHER = 0.01

smoothing_function(x) = (sqrt(x^2 + SMOOTHER^2) + x) / 2

data_points = collect(1:Δt:nrow(df))

function gr4j(params, output_times)

    X0 = params[1:nres+2]
    ODEparams = params[nres+3:end]

    time_span = (output_times[1], output_times[end])

    function odesystem!(dX, X, p, t)
        x1, x2, x3, x4 = p
        S, R = X[1], X[end]
        Sh = X[2:end-1]

        P = prec(t)
        Ep = pet(t)

        # Pn = max(P - Ep, 0.0)
        Pn = smoothing_function(P-Ep)
        Ps = smoothing_function(Pn * (1 - (S / x1)^α))
        # En = max(Ep - P, 0.0)
        En = smoothing_function(Ep - P)
        Es = smoothing_function(En * (2 * S / x1 - (S / x1)^α))
        Perc = x1^(1-β) / (Uₜ*(β-1)) * ν^(β-1) * S^β
        dS = Ps - Es - Perc

        Pr = Pn - Ps + Perc
        Qsh = (nres-1)*x4 * Sh[1:end]
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = x2 * (smoothing_function(R) / smoothing_function(x3))^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSh1
        dX[3:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt,TRBDF2())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

ODEparams = [1000.0, 2.0, 200.0, 2.5]
params = vcat(ones(nres+2) .* 50.0, ODEparams)
lower_bound = vcat(ones(nres+2), [1.0, -20.0, 1.0, 0.5])
upper_bound = vcat(ones(nres+2) .* 2000.0, [2000.0, 20.0, 300.0, 15.0])
times = 1.0:Δt:train_
Q = gr4j(params, times)

# UDE_model(params, output_times) = gr4jsnowNN(params, output_times, ann)

function NSE_loss(model, params, target_data, target_time)

    predicted_data = model(params, target_time)
    NSE(obs, pred) = 1 - sum((obs .- pred).^2) / sum((obs .- mean(obs)).^2)
    loss = -NSE(target_data, predicted_data)

    return loss 

end

function MSE_loss(model, params, target_data, target_time)

    predicted_data = model(params, target_time)
    MSE(obs, pred) = sum((obs .- pred).^2)
    loss = MSE(target_data, predicted_data)

    return loss

end

loss_function(p) = NSE_loss(gr4j, p, train_Y, times)
loss_function(params)

current_param = []

function callback(p, l)

    push!(current_param, p)
    println("NSE: "*string(-l))
    return false

end

opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoForwardDiff())

opt_problem = Optimization.OptimizationProblem(opt_func, params, lb=lower_bound, ub=upper_bound)
# opt_problem = Optimization.OptimizationProblem(opt_func, params)

# optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()
# optimizer = ADAM(0.1)
optimizer = CMAEvolutionStrategyOpt()
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=1000)

function gr4jsnow(params, output_times)

    X0 = params[1:nres+3]
    ODEparams = params[nres+4:end]

    time_span = (output_times[1], output_times[end])

    function odesystem!(dX, X, p, t)
        x1, x2, x3, x4, TRS, TRANS, DDF, Tbase = p
        S, R = X[1], X[end]
        Snow = X[2]
        Sh = X[3:end-1]

        P = prec(t)
        Ep = pet(t)
        T = temp(t)

        snowfrac = min(max((TRS + TRANS - T)/(2*TRANS), zero(T)), one(T))
        Psnow = P*snowfrac
        snowmelt = max(zero(Snow), min(DDF*(T-Tbase), Snow))
        dSnow = Psnow - snowmelt

        P = P - Psnow
        Pn = max(P - Ep, zero(P))
        Ps = Pn * (1 - (S / x1)^α)
        En = max(Ep - P, zero(P))
        Es = En * (2 * S / x1 - (S / x1)^α)
        Perc = x1^(1-β) / (Uₜ*(β-1)) * ν^(β-1) * S^β
        dS = Ps - Es - Perc + snowmelt

        Pr = Pn - Ps + Perc
        Qsh = (nres-1)/x4 * Sh[1:end]
        Quh = (nres-1)/x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * R^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSnow
        dX[3] = dSh1
        dX[4:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt,Tsit5())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q 
    
end

ODEparams = [1000.0, 2.0, 200.0, 2.5, 10.0, 13.0, 15.0, 0.0]
params = vcat(ones(nres+3), ODEparams)
lower_bound = vcat(ones(nres+3), [1.0, -20.0, 1.0, 0.5, -35.0, -30.0, 0.5, -30.0])
upper_bound = vcat(ones(nres+3) .* 2000.0, [2000.0, 20.0, 300.0, 15.0, 30.0, 35.0, 30.0, 30.0])

times = 1.0:Δt:train_
Qsnow, sol = gr4jsnow(params, times)

loss_function(p) = NSE_loss(gr4jsnow, p, train_Y, times)
loss_function(params)

opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoForwardDiff())

opt_problem = Optimization.OptimizationProblem(opt_func, params, lb=lower_bound, ub=upper_bound)
# opt_problem = Optimization.OptimizationProblem(opt_func, params)

# optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited()
# optimizer = ADAM(0.1)
optimizer = CMAEvolutionStrategyOpt()
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=1000)

#plot(times, sol[2,:])
#plot!(times, Qsnow./2)
#plot!(times,Q./2)


function gr4jsnowNN(params, output_times, ann)

    # X0 = params[1:nres+3]
    # X0 = params[1]
    X0 = params.ODE_states
    # ODEparams = params[nres+4:nres+7]
    # ODEparams = params[2]
    ODEparams = params.ODEparams
    # NNparams = params[end]

    time_span = (output_times[1], output_times[end])

    function odesystem!(dX, X, p, t)
        # NNparams = ComponentArray(p[end])
        NNparams = p.NNparams
        # ODEparams = p[nres+4:end-1]
        ODEparams = p.ODEparams
        x1, x2, x3, x4 = ODEparams
        S, R = X[1], X[end]
        Snow = X[2]
        Sh = X[3:end-1]

        P = prec(t)
        Ep = pet(t)
        T = temp(t)

        ann_outputs = ann([Snow, P, T], NNparams)
        snowfrac = ann_outputs[1]
        Psnow = P*snowfrac
        snowmelt = ann_outputs[2]
        dSnow = Psnow - snowmelt

        P = P - Psnow
        Pn = max(P - Ep, zero(P))
        Ps = Pn * (1 - (S / x1)^α)
        En = max(Ep - P, zero(P))
        Es = En * (2 * S / x1 - (S / x1)^α)
        Perc = x1^(1-β) / (Uₜ*(β-1)) * ν^(β-1) * S^β
        dS = Ps - Es - Perc + snowmelt

        Pr = Pn - Ps + Perc
        Qsh = (nres-1)/x4 * Sh[1:end]
        Quh = (nres-1)/x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * R^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSnow
        dX[3] = dSh1
        dX[4:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, params)
    sol = solve(prob, saveat = Δt,Tsit5())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

function initializeNN()
    rng = Random.default_rng()
    NNmodel = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16,2,tanh))
    p_NN_init, st_NN_init = Lux.setup(rng, NNmodel)
    ann(x, p) = NNmodel(x, p, st_NN_init)[1]
    # p_NN_init = ComponentArray(p_NN_init)

    return ann, p_NN_init

end

ann, initial_NN_params = initializeNN()
ODEparams = [1000.0, 2.0, 200.0, 2.5]
ODE_states = ones(nres+3)
# initial_params = (ODE_states, ODEparams, initial_NN_params)
initial_params = ComponentArray(ODE_states=ODE_states, ODEparams=ODEparams, NNparams=initial_NN_params)

gr4jsnowNN(initial_params, train_points, ann)
UDE_model(params, output_times) = gr4jsnowNN(params, output_times, ann)

function NSE_loss(model, params, target_data, target_time)

    predicted_data = model(params, target_time)
    NSE(obs, pred) = 1 - sum((obs .- pred).^2) / sum((obs .- mean(obs)).^2)
    loss = -NSE(target_data, predicted_data)

    return loss 

end

loss_function(p) = NSE_loss(UDE_model, p, train_Y, train_points)
loss_function(initial_params)

function callback(p, l)

    println("NSE: "*string(-l))
    return false

end

opt_func = Optimization.OptimizationFunction((p, known_params) -> loss_function(p), Optimization.AutoZygote())

opt_problem = Optimization.OptimizationProblem(opt_func, initial_params)

optimizer = ADAM(0.1)
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=10)

out_params = sol.u
times = 1.0:Δt:train_

Q_nn = gr4jsnowNN(out_params, times, ann)

plot(times, Q_nn)
plot!(times, df[1:train_,:streamflow], dpi = 300)
plot!(times,Q_nn)

#test data
times = train_ + 1:Δt:data_points[end]

Q_nn = gr4jsnowNN(out_params, times, ann)

plot(times, Q_nn)
plot!(times, df[(train_ +1):end,:streamflow], dpi = 300)
plot!(times,Q_nn)




savefig("chepe_plot.png")
