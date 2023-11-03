#=
    This module describes the EC coupling starting from the framework of the
    Shannon-Bers model, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#
using DifferentialEquations
using ModelingToolkit

export get_ecc_equations

function beta_cai(Cai)
    @parameters TrpnTotal = 35
    @parameters KmTrpn = 0.5
    @parameters CmdnTotal = 50
    @parameters KmCmdn = 2.38
    return inv(1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)
end

function get_ecc_equations()
    @variables t
    D = Differential(t)

    @variables i_CaJSR(t) CaNSR(t) V(t) Nai(t) Ki(t) i_b(t) i_g(t) i_d(t) i_f(t) i_fca(t) i_y(t) i_r(t) i_s(t) i_sslow(t)
    @variables i_Nam(t) i_Nah(t) i_Naj(t) i_PO1(t) i_PO2(t) i_PC2(t) i_nKs(t) i_CK1(t) i_CK2(t) i_OK(t) i_IK(t)
    @variables Cai_sub_SR(t) Cai_sub_SL(t)

    # -------------------------------------------------------------------------
    # Model of ECC of rat neonatal ventricular myocyte 2009
    # Code & model: Topi Korhonen, University of Oulu (topi.korhonen@oulu.fi)
    # 
    # PLEASE MENTION THE FOLLOWING REFERENCE WHEN USING THIS CODE OR PART OF IT:
    # Korhonen et al. "Model of excitation-contraction coupling of rat neonatal 
    # ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
    #
    # ONLY FOR ACADEMIC USE, DO NOT DISTRIBUTE
    # -------------------------------------------------------------------------
    
    # Index numbers, put the Cai equations to the end -> increase i_Cai_sub_SR when
    # adding other odes
    
    
    # Physical constants
    F 	= 96.5
    T 	= 305
    R 	= 8.314
    Cm 	= 1.0
    
    # Ion concentrations in DMEM
    Cao = 1796
    Nao = 154578
    Ko = 5366
    
    # Cell geometry
    rSR_true = 6
    rSL_true = 10.5
    
    # Ca diffusion grid
    dx = 0.1
    rSR = rSR_true + 0.5*dx
    rSL = rSL_true - 0.5*dx
    j = round(rSR/dx):1:round(rSL/dx) # Spatial index of Cai diffusion
    #j = j'
    
    # More cell geometry
    V_sub_SR = 4/3 * pi * (rSR_true + dx)^3 / 1000 - 4/3 * pi * (rSR_true)^3 / 1000 # pl
    V_sub_SL = 4/3 * pi * rSL_true^3 / 1000 - 4/3 * pi * (rSL_true-dx)^3 / 1000 #pl
    Acap = 4 * pi * rSL_true^2 * 1e-8 # cm^2
    VSR = 0.043*1.5*1.4
    VNSR = 0.9*VSR
    VJSR = VSR-VNSR
    Vmyo = 4/3 * pi * rSL_true^3 / 1000 - 4/3 * pi * rSR_true^3 / 1000
    

    # Ca buffers
    csqntot = 24750
    Kmcsqn = 800
    betaSR = inv( 1 + csqntot*Kmcsqn./(i_CaJSR + Kmcsqn).^2)
    
    # Diffusion coefficient
    D = 7 # mum^2/ms set to achive correct diff. speed 0.31 mum/ms
    
    # NCX from Pandit rat model
    fNaCa = 1
    kNaCa = 2.2680e-016
    dNaCa = 1e-16
    gamma = 0.5
    
    # INaK 
    INaKmax = 2.7
    KmNai 	= 18600
    nNaK    = 3.2
    KmKo 	= 1500
    
    
    # --------------------------------------------------------
    # Equations used in odes 
    # --------------------------------------------------------
    
    ENa = R*T/F * log( (Nao) / (Nai) )
    EK = R*T/F * log( (Ko) / (Ki) )
    ECa = R*T/2/F * log( Cao / Cai_sub_SL)
    
    # NCX
    INaCa = kNaCa * ( (exp( 0.03743*gamma.*V ) .* Nai.^3 .* Cao - exp( 0.03743*(gamma-1).*V ) .* Nao^3 .* Cai_sub_SL .* fNaCa ) / ( 1 + dNaCa*(Nao^3 .* Cai_sub_SL.*fNaCa + Nai.^3 .* Cao) ) )
    
    # L-type calcium current 
    GCaL 	= 1.3125e-4*0.8*0.6
    ICaL = GCaL * i_d * i_f * i_fca * 4 * V * F^2 / R / T * ( Cai_sub_SL * exp(2*V*F/R/T) 
        - 0.341 * Cao) / ( exp(2*V*F/R/T) - 1)
    dinf = 1 / ( 1 + exp( (11.1+V)/-7.2 ) )
    alphad = 1.4 / ( 1 + exp( (-35-V)/13 ) ) + 0.25
    betad = 1.4 / ( 1 + exp( (V + 5)/5 ) )
    gammad = 1 / ( 1 + exp( (50 - V)/20 ) )
    taud = alphad * betad + gammad
    finf = 1 / ( 1 + exp( (V + 23.3)/5.4 ) )
    tauf = 1125 * exp( -(V + 27)^2/240 ) + 165/(1 + exp( (25-V)/10 )) + 120
    alphafca = 1 / ( 1 + (Cai_sub_SL*1e-3/(0.000325*1.5))^8 )
    betafca = 0.1 / ( 1 + exp( (Cai_sub_SL*1e-3 - 0.0005)/0.0001 ) )
    gammafca = 0.2 / ( 1 + exp( (Cai_sub_SL*1e-3 - 0.00075)/0.0008 ) )
    fcainf = (alphafca + betafca + gammafca + 0.23)/1.46
    taufca = 10 # modif 
    kfca = ifelse( ( fcainf > i_fca ) & ( V > (-60) ), 0, 1)
    
    # T-Type
    gCaT = 0.2
    binf = 1 / ( 1 + exp( -( V + 37.49098) / 5.40634 ) )
    taub = 0.6 + 5.4 / ( 1 + exp( ( V + 100 ) * 0.03 ) )
    ginf = 1 / ( 1 + exp( ( V + 66 ) / 6 ) )
    taug = 1 + 40 / ( 1 + exp( ( V + 65 ) * 0.08 ) )
    ICaT = gCaT * i_b .* i_g * ( V - ECa + 106.5)
    
    # Cab
    gCab = 0.0008
    ICab = gCab * ( V - ECa)
    
    # Nab
    gNab = 0.0026
    INab = gNab * ( V - ENa)
    
    # If 
    gf = 0.021
    fNa = 0.2
    fK = 1-fNa
    yinf = 1 ./ (1 + exp( (V+78.65)./6.33 ) ) # Fitted
    tauy = 1 ./ ( 0.11885.*exp( (V+75)./28.37 ) + 0.56236.*exp( (V+75) ./ -14.19 ) ) .* 1000
    IfNa = gf*i_y* fNa*(V-ENa)
    IfK = gf*i_y*fK*(V-EK)
    If = gf*i_y*( fNa*(V-ENa) + fK*(V-EK) )
    
    # IK1
    IK1 = 0.0515 .* ( Ko/(Ko + 210) ) .* ( ( V - EK -6.1373 ) ./ ( 0.1653  + exp( 0.0319 * ( V - EK - 6.1373) ) ) )
    
    # Ito
    gt = 0.1
    sinf = 1 ./ (1 + exp( (V+31.97156)./4.64291 ) ) 
    rinf = 1 ./ (1 + exp( (V-3.55716)./-14.61299 ) ) 
    slowinf = sinf
    taur = 1 ./ ( 45.16 .* exp( 0.03577.*(V + 50)) + 98.9 .* exp(-0.1.*(V + 38)) ) * 1000
    taus = ( 0.35 .* exp( -( ((V+70)/15).^2) ) + 0.035 ) * 1000-26.9
    tausslow = ( 3.7 .* exp( -( ((V+70)/30).^2) ) + 0.035 ) * 1000+37.4
    Ito = gt * i_r .* (0.706 .* i_s + 0.294 .* i_sslow ) * ( V - EK )
    
    # INa
    gNa = 35
    Naminf = 1/(1+exp((V+45)/-6.5))
    Nahinf = 1/(1+exp((V+76.1)/6.07))
    Najinf = Nahinf
    Nataum = 0.00136 / ( 0.32*(V + 47.13)/(1 - exp(-0.1*(V + 47.13))) + 0.08*exp(-V/11) )
    Natauh = ifelse(V >= -40, 0.0004537 .* (1 + exp( (V+10.66)./-11.1 )), 0.00349 ./ ( 0.135.*exp( (V+80)./-6.8) + 3.56.*exp(0.079.*V) + 3.1e5.*exp(0.35.*V)))
    Natauj = ifelse(V >= -40, 0.01163 .* (1 + exp( -0.1.*(V+32))) ./ exp(-2.535e-7.*V), 0.00349 ./ ( (V + 37.78)./(1 + exp(0.311.*(V + 79.23))) .* (-127140 .*exp(0.2444.*V) - 3.474e-5.*exp(-0.04391.*V) ) + 0.1212.*exp(-0.01052.*V)./(1 + exp(-0.1378.*(V + 40.14))) ))
    
    INa = gNa*i_Nam.^3*i_Nah*i_Naj*( V - ENa )

    # RyR
    n = 4
    PC1 = 1 - i_PO1
    nu1 = 0.01*3/3
    KmRyR =1.35*2.6./(1+exp((i_CaJSR-530)./200))+1.5-0.9-0.3-0.05
    kapos = 3/3
    kaneg = 0.48/3
    Jrel = nu1 * ( i_PO1 ) * ( i_CaJSR - Cai_sub_SR )
    
    # SR Ca-ATPase 
    Vmaxf = 0.9996
    Vmaxr = Vmaxf
    Kmf = 0.5
    Kmr = 7000*Kmf
    Hf = 2
    Hr =1*Hf
    k = 5e-6
    
    
    Jup = ( Vmaxf.*(Cai_sub_SR./Kmf).^Hf - Vmaxr.*(CaNSR./Kmr).^Hr ) ./ ( 1 + (Cai_sub_SR./Kmf).^Hf + (CaNSR./Kmr).^Hr )
    
    # Jleak
    Jleak = k*(CaNSR-Cai_sub_SR)
    
    # IKs
    GKs 	= 0.05
    IKs = GKs * i_nKs^2 * ( V - EK )
    alphan = 0.00000481333 * ( V + 26.5 ) / ( 1 - exp( -0.128*( V + 26.5 ) ) )
    betan = 0.0000953333 * exp( -0.038 * ( V + 26.5 ) )
    nKsinf = alphan / ( alphan + betan )
    nKstau = 750
    
    # Na/K pump current
    sigma = 1/7 * ( exp( Nao/67300 ) - 1)
    fNaK = 1 / ( 1 + 0.1245 * exp( -0.1*V*F/R/T ) + 0.0365 * sigma * exp( -V*F/R/T ) )
    INaK = INaKmax * fNaK / ( 1 + ( KmNai/Nai )^(nNaK) ) * Ko / ( Ko + KmKo)
    
    # IKr
    GKr 	= 0.06
    kf  = 0.023761
    kb 	= 0.036778
    IKr = i_OK * GKr * ( V - R*T/F * log( ( 0.98*Ko + 0.02*Nao ) / ( 0.98*Ki + 0.02*Nai ) ) )
    CK0 = 1 - ( i_CK1 + i_CK2 + i_OK + i_IK )
    alphaa0 = 0.022348 * exp( 0.01176*V )
    betaa0 = 0.047002 * exp( -0.0631*V )
    alphaa1 = 0.013733 * exp( 0.038198*V )
    betaa1 = 0.0000689 * exp( -0.04178*V )
    alphai_mERG = 0.090821 * exp( 0.023391 * V )
    betai_mERG = 0.006497 * exp( -0.03268 * V )
    
    Jtr = (CaNSR - i_CaJSR)/ 200
    
    # Ca fluxes
    JCa_SL = (2 * INaCa - ICaL - ICaT - ICab) * Acap*Cm/2/1/F*1e6
    JCa_SR = Jleak - Jup + Jrel
    
    # -------------------------------------------------------------------
    # Differential equations 
    # -------------------------------------------------------------------
    # Cai
    m = length(j)
    @variables Cai(t)[1:m]
    @parameters Dca = 7 # mum^2/ms set to achive correct diff. speed 0.31 mum/ms

    #=
    @parameters TrpnTotal = 35
    @parameters KmTrpn = 0.5
    @parameters CmdnTotal = 50
    @parameters KmCmdn = 2.38
    #beta_cai = (1 + TrpnTotal * KmTrpn / (Cai + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai + KmCmdn)^2)^(-1)
    for i in 1:m
        beta_eq = [
           D(beta_cai[i]) ~ inv(1 + TrpnTotal * KmTrpn / (Cai[i] + KmTrpn)^2 + CmdnTotal * KmCmdn / (Cai[i] + KmCmdn)^2)
        ]
        push!(beta_eqs, beta_eq)
    end
    =#

    eq = [
        D(Cai[1]) ~ (Dca / (j[1] * dx^2) * ((1 + j[1]) * Cai[2] - 2 * j[1] * Cai[1] + (j[1] - 1) * Cai[1] ) + JCa_SR/V_sub_SR) * beta_cai(Cai[1])
    ]

    eqs = [eq, Cai_sub_SR~Cai[1], Cai_sub_SL~Cai[m]]

    for i in 2:m-1
        eq = [
           D(Cai[i]) ~ (Dca / (j[i] * dx^2) * ((1 + j[i]) * Cai[i+1] - 2 * j[i] * Cai[i] + (j[i] - 1) * Cai[i-1] )) * beta_cai(Cai[i])
        ]
        push!(eqs, eq)
    end

    eq = [
        D(Cai[m]) ~ (Dca / (j[m] * dx^2) * ((1 + j[m]) * Cai[m] - 2 * j[m] * Cai[m] + (j[m] - 1) * Cai[m-1] ) + JCa_SL//V_sub_SL) * beta_cai(Cai[m])
    ]
    push!(eqs, eq)


    # Istim
    Istim = ifelse(mod(t,cycleLength) <= 2, -80, 0.0)
    
    # Other odes
    SR_eqs = [
        D(i_CaJSR) ~ betaSR * (-Jrel + Jtr)/VJSR,
        D(CaNSR) ~ ( Jup - Jleak - Jtr)/VNSR,
        D(V) ~ -( INab + INaCa + ICaL + ICaT + If + Ito + IK1 + IKs + IKr + INa + INaK + ICab +Istim ) / Cm
    ]
    # Nai, Ki
    Na_K_eqs = [
        D(Nai) ~ -( IfNa + INab + INa + 3 * INaCa + 3 * INaK ) * Acap*Cm/F*1e6/Vmyo,
        D(Ki) ~ -( IfK + Ito + IK1 + IKs + IKr + Istim  - 2 * INaK ) * Acap*Cm/F*1e6/Vmyo
    ]
    
    # T-Type
    Ttype_eqs = [
        D(i_b) ~ (binf - i_b ) / taub,
        D(i_g) ~ (ginf - i_g ) / taug
    ]

    # L-type calcium current
    Ltype_eqs = [
        D(i_d) ~ ( dinf - i_d ) / taud,
        D(i_f) ~ ( finf - i_f ) / tauf,
        D(i_fca) ~ kfca * ( fcainf - i_fca ) / taufca
    ]

    
    # If
    IF_eqs = [
        D(i_y) ~ (yinf - i_y)/tauy
    ]
    
    
    # IKto
    IKto_eqs = [
        D(i_r) ~ (rinf - i_r ) / taur,
        D(i_s) ~ (sinf - i_s ) / taus,
        D(i_sslow) ~ (slowinf - i_sslow ) / tausslow
    ]
    
    # INa
    INa_eqs = [
        D(i_Nam) ~ (Naminf - i_Nam)/Nataum/1000,
        D(i_Nah) ~ (Nahinf - i_Nah)/Natauh/1000,
        D(i_Naj) ~ (Najinf - i_Naj)/Natauj/1000
    ]

    
    # RyR
    ryr_eqs = [
        D(i_PO1) ~ kapos * Cai_sub_SR^n / (Cai_sub_SR^n + KmRyR^n) * PC1 - kaneg * i_PO1,
        D(i_PO2) ~ 0,
        D(i_PC2) ~ 0
    ]

    
    # IKs
    IKs_eqs = [
       D(i_nKs) ~ ( nKsinf - i_nKs ) / nKstau 
    ]
    
    
    # Rapid delayed rectifier K current (mERG)
    IKr_eqs = [
        D(i_CK1) ~ alphaa0 * CK0 - betaa0 * i_CK1 + kb * i_CK2 - kf * i_CK1,
        D(i_CK2) ~ kf * i_CK1 - kb * i_CK2 + betaa1 * i_OK - alphaa1 * i_CK2,
        D(i_OK) ~ alphaa1 * i_CK2 - betaa1 * i_OK + betai_mERG * i_IK - alphai_mERG * i_OK,
        D(i_IK) ~ alphai_mERG * i_OK - betai_mERG * i_IK
    ]


    return vcat(eqs, SR_eqs, Na_K_eqs, Ttype_eqs, Ltype_eqs, IF_eqs, IKto_eqs, INa_eqs, INa_eqs, ryr_eqs, IKs_eqs, IKr_eqs)
end


