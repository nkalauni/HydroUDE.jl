# Load all external packages

using Revise

using DataFrames, Dates, Statistics
using CSV

using DifferentialEquations
# using OrdinaryDiffEq
using Lux, DiffEqFlux
using ComponentArrays
using SciMLSensitivity

using Optimization, OptimizationBBO
using OptimizationOptimisers
using Zygote

using Interpolations

using Plots

using Random
Random.seed!(123)

df = CSV.read("src/chepe_data.csv", DataFrame; dateformat="mm/dd/yyyy")

interpolate_method = SteffenMonotonicInterpolation()

const Δt = 1.0
data_points = collect(1:nrow(df))

temp = interpolate(data_points, df[:,:temperature], interpolate_method)
prec = interpolate(data_points, df[:,:precipitation], interpolate_method)
pet = interpolate(data_points, df[:,:pet], interpolate_method)
flow = df[:, :streamflow]

const split_ratio = 0.7
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

        Pn = max(P - Ep, zero(P))
        Ps = Pn * (1 - (S / x1)^α)
        En = max(Ep - P, zero(P))
        Es = En * (2 * S / x1 - (S / x1)^α)
        Perc = x1^(1-β) / (Uₜ*(β-1)) * ν^(β-1) * S^β
        dS = Ps - Es - Perc

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
        dX[2] = dSh1
        dX[3:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt)

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

ODEparams = [1000.0, 2.0, 200.0, 2.5]
params = vcat(ones(nres+2), ODEparams)
times = 1.0:Δt:1095.0
Q = gr4j(params, times)

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
    sol = solve(prob, saveat = Δt)

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q, sol
    
end

ODEparams = [1000.0, 2.0, 200.0, 2.5, 10.0, 13.0, 15.0, 0.0]
params = vcat(ones(nres+3), ODEparams)
times = 1.0:Δt:1095.0
Qsnow, sol = gr4jsnow(params, times)

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
    sol = solve(prob, saveat = Δt)

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

optimizer = ADAM(0.001)
sol = Optimization.solve(opt_problem, optimizer, callback=callback, maxiters=50)

