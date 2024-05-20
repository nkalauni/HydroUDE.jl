# Define model constants for GR4J model
const α = 2
const β = 5
const γ = 5
const ω = 3.5
const ϵ = 1.5
const Φ = 0.9
const ν = 4/9
const Uₜ = 1
const nres = 11

function gr4j(params, output_times, args...)

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
        Qsh = (nres-1)*x4 * Sh[1:end]   # x4 = 1/x4
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * max(R,0)^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSh1
        dX[3:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt, Tsit5())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

function gr4jsnow(parameters, output_times, args...)

    if length(args)>0
        X0 = args[1]    #initial_parameters
        ODEparams = parameters
    else
        X0 = parameters[1:nres+3]
        ODEparams = parameters[nres+4:end]
    end

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
        Qsh = (nres-1)*x4 * Sh[1:end]   #x4 = 1/x4  --------
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * max(0,R)^ω
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


function GR4J_model(parameters, output_times, args...)

    if length(args)>0
        X0 = args[1]    #initial_parameters
        ODEparams = parameters
    else
        X0 = parameters[1:nres+2]
        ODEparams = parameters[nres+3:end]
    end
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
        Qsh = (nres-1)*x4 * Sh[1:end]   # x4 = 1/x4
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * max(R,0)^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSh1
        dX[3:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt, Rosenbrock23())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

function GR4J_model()
    return "GR4J_model"
end

"""
Function definiting gr4j model with NN for snow component.
"""
function initializeNN()
    rng = Random.default_rng()
    NNmodel = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16,2,tanh))
    p_NN_init, st_NN_init = Lux.setup(rng, NNmodel)
    ann(x, p) = NNmodel(x, p, st_NN_init)[1]
    # p_NN_init = ComponentArray(p_NN_init)

    return ann, p_NN_init

end


# --

function gr4jsnowNN(parameters, output_times, ann, args...)
    
    # if length(args)>0
    #     ann = args[1]
    # else
    #     print("यहाँले न्युरल नेटवर्क को मोडेल पठाउन बिर्सिनिनु भयो जस्तो छ, कृपया आर्गुमेंट पुन: निरिक्षण गर्नुहोस। ")
    # end


    if length(args)>0
        X0 = args[1]    #initial_parameters
    else
        X0 = parameters.ODEstates
    end
    ODEparams = parameters.ODEparams
    NNparams = parameters.NNparams

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
        Qsh = (nres-1)*x4 * Sh[1:end]   #x4 = 1/x4
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * max(0,R)^ω    # avoid complex number
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSnow
        dX[3] = dSh1
        dX[4:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, parameters)
    sol = solve(prob, saveat = Δt,Tsit5())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end

function gr4jSG(parameters, output_times, args...)

    if length(args)>0
        X0 = args[1]    #initial_parameters
        ODEparams = parameters
    else
        X0 = parameters[1:nres+4]
        ODEparams = parameters[nres+5:end]
    end

    time_span = (output_times[1], output_times[end])

    function odesystem!(dX, X, p, t)
        x1, x2, x3, x4, TRS, TRANS, DDF_snow, Tbase, DDF_ice = p
        S, R = X[1], X[end]  
        Snow = X[2]
        Ice = X[3]
        Sh = X[4:end-1]

        P = prec(t)
        Ep = pet(t)
        T = temp(t)

        snowfrac = min(max((TRS + TRANS - T)/(2*TRANS), zero(T)), one(T))
        Psnow = P*snowfrac
        snowmelt = max(zero(Snow), min(DDF_snow*(T-Tbase), Snow))
        dSnow = Psnow - snowmelt

        #TODO icemelt only if snow is zero ---> how to implement this?
        icemelt = max(zero(Ice), min(DDF_ice*(T-Tbase), Ice))
        dIce = -icemelt

        P = P - Psnow
        Pn = max(P - Ep, zero(P))
        Ps = Pn * (1 - (S / x1)^α)
        En = max(Ep - P, zero(P))
        Es = En * (2 * S / x1 - (S / x1)^α)
        Perc = x1^(1-β) / (Uₜ*(β-1)) * ν^(β-1) * S^β
        dS = Ps - Es - Perc + snowmelt

        Pr = Pn - Ps + Perc
        Qsh = (nres-1)*x4 * Sh[1:end]   #x4 = 1/x4  --------
        Quh = (nres-1)*x4 * Sh[end]
        dSh1 = Pr - Qsh[1]
        dSh = [Qsh[i] - Qsh[i+1] for i in 1:nres-1]

        Q9 = Φ*Quh
        F = (x2 / x3^ω) * max(0,R)^ω
        Qr = x3^(1-γ)/(Uₜ*(γ-1)) * R^γ
        dR = Q9 + F - Qr

        dX[1] = dS
        dX[2] = dSnow
        dX[3] = dIce
        dX[4] = dSh1
        dX[5:end-1] = dSh
        dX[end] = dR

    end

    prob = ODEProblem(odesystem!, X0, time_span, ODEparams)
    sol = solve(prob, saveat = Δt, Tsit5())

    Quh = (nres-1)/ODEparams[4] .* sol[nres+1,:]
    F = (ODEparams[2] / ODEparams[3]^ω) .* sol[end,:]
    Qr = ODEparams[3]^(1-γ) / (Uₜ*(γ-1)) .* sol[end,:].^γ
    Qd = max.(0, (1-Φ) .* Quh .- F)
    Q = Qr + Qd

    return Q
    
end