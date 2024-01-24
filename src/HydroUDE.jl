# Load all external packages

using Revise

using DataFrames, Dates, Statistics
using CSV

using DifferentialEquations
using Lux, DiffEqFlux
using SciMLSensitivity

using Optimization, OptimizationBBO
using Zygote

using Interpolations

using Plots

df = CSV.read("src/chepe_data.csv", DataFrame; dateformat="mm/dd/yyyy")
data_points = collect(1:nrow(df))

interpolate_method = SteffenMonotonicInterpolation()

temp = interpolate(data_points, df[:,:temperature], interpolate_method)
prec = interpolate(data_points, df[:,:precipitation], interpolate_method)
pet = interpolate(data_points, df[:,:pet], interpolate_method)

function full_hydrograph(basetime; timestep = 0.05)
    times = collect(1:timestep:2*basetime)
    Uhydrograph = zeros(length(times))
    Scurve = zeros(length(times)+1)

    for (i, t) in enumerate(times)
        if t <= basetime
            Scurve[i+1] = 0.5 * (t / basetime)^2.5
        elseif basetime < t <= 2*basetime
            Scurve[i+1] = 1.0 - 0.5 * (2.0 - t / basetime)^2.5
        end
    end

    Uhydrograph = Scurve[2:end] - Scurve[1:end-1]
    return Uhydrograph

end

function half_hydrograph(basetime; timestep = 0.05)
    times = collect(1:timestep:basetime)
    Uhydrograph = zeros(length(times))
    Scurve = zeros(length(times)+1)

    for (i, t) in enumerate(times)
        if t <= basetime
            Scurve[i+1] = 0.5 * (t / basetime)^2.5
        elseif basetime < t <= 2*basetime
            Scurve[i+1] = 1.0 - 0.5 * (2.0 - t / basetime)^2.5
        end
    end

    Uhydrograph = Scurve[2:end] - Scurve[1:end-1]
    return Uhydrograph

end

function gr4j(params, output_times)

    function ddesystem!(dS, S, h, p, t)
        x1, x2, x3, x4 = p
        soil_store, routing_store = S

        P = prec(t)
        Ep = pet(t)

        Pn = max(P - Ep, zero(P))
        Ps = Pn * (1 - (soil_store / x1)^2)
        En = max(Ep - P, zero(P))
        Es = En * (2 * Soil_store / x1 - (soil_store / x1)^2)
        Perc = x1^(-4) / 4 * (4/9)^4 * soil_store^5
        dS[1] = Ps - Es - Perc

        ground_flow = 0.9 * (Pn - Ps + Perc)
        direct_runoff = 0.1 * (Pn - Ps + Perc)
        Q90 = route(ground_flow, half_uh)
        Q10 = route(direct_runoff, full_uh)
        recharge = x2 * ((max(routing_store, zero(routing_store)) / x3)^3.5)
        Qrouted = x3^(-4) / 4 * routing_store^5
        dS[2] = Q90 + recharge - Qrouted

        prob = DDEProblem

    end

end


const alpha = 2; const beta = 5; 

const u_0 = 0.0
const t_0 = 0.0
const T = 5.0
const n = 1000
const delta_t = 0.05

const lags = range(delta_t, length=n, stop=(T-t_0)*(1-1/n))

# left Riemann sum
function f(du, u, h, p, t)
    alpha, beta = p
    du[1] = 1 - alpha*u[1] - beta*((T-t_0)/n)*sum([h(p, t-tau)[1] for tau in lags])
end

h(p, t) = zeros(1)
tspan = (0.0, T)
u0 = [u_0]
prob = DDEProblem(f, u0, h, tspan, [alpha, beta])
algo = MethodOfSteps(Tsit5())
sol = solve(prob, algo)
plot(sol)