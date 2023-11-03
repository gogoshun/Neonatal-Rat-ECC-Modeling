#=
    This module describes the beta-adrenergic signaling pathway in mouse
    ventricular myocyte, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#
using ModelingToolkit

export get_bar_equations

function get_bar_equations()
    @variables t
    D = Differential(t)

    @variables LR(t) LRG(t) RG(t) b1AR_S464(t) b1AR_S301(t) GsaGTPtot(t) GsaGDP(t) Gsby(t) AC_GsaGTP(t) PDEp(t) 
    @variables cAMPtot(t) RC_I(t) RCcAMP_I(t) RCcAMPcAMP_I(t) RcAMPcAMP_I(t) PKACI(t) PKACI_PKI(t) RC_II(t) RCcAMP_II(t) RCcAMPcAMP_II(t)
    @variables RcAMPcAMP_II(t) PKACII(t) PKACII_PKI(t) I1p_PP1(t) I1ptot(t) PLBp(t) PLMp(t) LCCap(t) LCCbp(t) RyRp(t)
    @variables TnIp(t) KS79(t) KS80(t) KSp(t) CFTRp(t) KURp(t)

    # Drug concentrations
    Ligtot = 0.0               # [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
    FSK=0
    IBMX=0
    LCCtotBA = 0.025           # [uM] - [umol/L cytosol]
    plb_val=106                # MOUSE
    PP1_PLBtot = 0.89          # [uM] - [umol/L cytosol]
    PLMtotBA = 48              # [uM] - [umol/L cytosol] MOUSE
    PLBtotBA = plb_val                     # [uM] - [umol/L cytosol]
    ISO=Ligtot

    ## b-AR module
    b1ARtot         = 0.00528        # (uM) total b1-AR protein # MOUSE
    #b1ARtot=0.028  # RABBIT
    kf_LR           = 1              # (1/[uM ms]) forward rate for ISO binding to b1AR
    kr_LR           = 0.285          # (1/ms) reverse rate for ISO binding to b1AR
    kf_LRG          = 1              # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
    kr_LRG          = 0.062          # (1/ms) reverse rate for ISO:b1AR association with Gs
    kf_RG           = 1              # (1/[uM ms]) forward rate for b1AR association with Gs
    kr_RG           = 33             # (1/ms) reverse rate for b1AR association with Gs
    Gstot           = 3.83           # (uM) total Gs protein
    k_G_act         = 16e-3          # (1/ms) rate constant for Gs activation
    k_G_hyd         = 0.8e-3         # (1/ms) rate constant for G-protein hydrolysis
    k_G_reassoc     = 1.21           # (1/[uM ms]) rate constant for G-protein reassociation
    kf_bARK         = 1.1e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
    kr_bARK         = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
    kf_PKA          = 3.6e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
    kr_PKA          = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by PKA

    b1ARact = b1ARtot - b1AR_S464 - b1AR_S301
    b1AR = b1ARact - LR - LRG - RG
    Gs = Gstot - LRG - RG - Gsby
    bARK_desens = kf_bARK*(LR+LRG)
    bARK_resens = kr_bARK*b1AR_S464
    PKA_desens = kf_PKA*PKACI*b1ARact
    PKA_resens = kr_PKA*b1AR_S301
    G_act = k_G_act*(RG+LRG)
    G_hyd = k_G_hyd*GsaGTPtot
    G_reassoc = k_G_reassoc*GsaGDP*Gsby

    bar_eqs = [
        D(LR) ~ (kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs),    # du[1]
        D(LRG) ~ (kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG),                 # du[2]
        D(RG) ~ (kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG),                    # du[3]
        D(b1AR_S464) ~ (bARK_desens - bARK_resens),                         # du[4]
        D(b1AR_S301) ~ (PKA_desens - PKA_resens),                           # du[5]
        D(GsaGTPtot) ~ (G_act - G_hyd),                                     # du[6]
        D(GsaGDP) ~ (G_hyd - G_reassoc),                                    # du[7]
        D(Gsby) ~ (G_act - G_reassoc)                                       # du[8]
    ]

    ## cAMP module
    ACtot           = 70.57e-3        # (uM) total adenylyl cyclase # MOUSE
    # ACtot=47e-3  # RABBIT
    ATP             = 5e3             # (uM) total ATP
    k_AC_basal      = 0.2e-3          # (1/ms) basal cAMP generation rate by AC
    Km_AC_basal     = 1.03e3          # (uM) basal AC affinity for ATP
    
    Kd_AC_Gsa       = 0.4             # (uM) Kd for AC association with Gsa
    kf_AC_Gsa       = 1               # (1/[uM ms]) forward rate for AC association with Gsa
    kr_AC_Gsa       = Kd_AC_Gsa       # (1/ms) reverse rate for AC association with Gsa
    
    k_AC_Gsa        = 8.5e-3          # (1/ms) basal cAMP generation rate by AC:Gsa
    Km_AC_Gsa       = 315.0           # (uM) AC:Gsa affinity for ATP
    
    Kd_AC_FSK       = 44.0            # (uM) Kd for FSK binding to AC
    k_AC_FSK        = 7.3e-3          # (1/ms) basal cAMP generation rate by AC:FSK
    Km_AC_FSK       = 860.0           # (uM) AC:FSK affinity for ATP
    
    PDEtot          = 22.85e-3        # (uM) total phosphodiesterase
    k_cAMP_PDE      = 5e-3            # (1/ms) cAMP hydrolysis rate by PDE
    k_cAMP_PDEp     = 2*k_cAMP_PDE    # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
    Km_PDE_cAMP     = 1.3             # (uM) PDE affinity for cAMP
    
    Kd_PDE_IBMX     = 30.0            # (uM) Kd_R2cAMP_C for IBMX binding to PDE
    k_PKA_PDE       = 7.5e-3          # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
    k_PP_PDE        = 1.5e-3          # (1/ms) rate constant for PDE dephosphorylation by phosphatases

    cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II)
    AC = ACtot-AC_GsaGTP
    GsaGTP = GsaGTPtot - AC_GsaGTP
    AC_FSK = FSK*AC/Kd_AC_FSK
    AC_ACT_BASAL = k_AC_basal*AC*ATP/(Km_AC_basal+ATP)
    AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP)
    AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP)
    PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX
    PDE = PDEtot - PDE_IBMX - PDEp
    PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP)

    cAMP_eqs = [
        D(AC_GsaGTP) ~ (kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP),
        D(PDEp) ~ (k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp),
        D(cAMPtot) ~ (AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT)
    ]


    ## PKA module

    PKItot          = 0.18             # (uM) total PKI
    kf_RC_cAMP      = 1                # (1/[uM ms]) Kd for PKA RC binding to cAMP
    kf_RCcAMP_cAMP  = 1                # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
    kf_RcAMPcAMP_C  = 4.375            # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
    kf_PKA_PKI      = 1                # (1/[uM ms]) Ki for PKA inhibition by PKI
    kr_RC_cAMP      = 1.64             # (1/ms) Kd for PKA RC binding to cAMP
    kr_RCcAMP_cAMP  = 9.14             # (1/ms) Kd for PKA RC:cAMP binding to cAMP
    kr_RcAMPcAMP_C  = 1                # (1/ms) Kd for PKA R:cAMPcAMP binding to C
    kr_PKA_PKI      = 2e-4             # (1/ms) Ki for PKA inhibition by PKI
    epsilon         = 10               # (-) AKAP-mediated scaling factor

    PKI = PKItot - PKACI_PKI - PKACII_PKI

    PKA_eqs = [
        D(RC_I) ~ (- kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I),
        D(RCcAMP_I) ~ (- kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I),
        D(RCcAMPcAMP_I) ~ (- kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI),
        D(RcAMPcAMP_I) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I),
        D(PKACI) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI),
        D(PKACI_PKI) ~ (- kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI),
        D(RC_II) ~ (- kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II),
        D(RCcAMP_II) ~ (- kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II),
        D(RCcAMPcAMP_II) ~ (- kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII),
        D(RcAMPcAMP_II) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II),
        D(PKACII) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI),
        D(PKACII_PKI) ~ (- kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI)
    ]


    ## I-1/PP1 module
    I1tot           = 0.3             # (uM) total inhibitor 1
    k_PKA_I1        = 60e-3           # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
    Km_PKA_I1       = 1.0             # (uM) Km for I-1 phosphorylation by type 1 PKA
    Vmax_PP2A_I1    = 14.0e-3         # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
    Km_PP2A_I1      = 1.0             # (uM) Km for I-1 dephosphorylation by PP2A
    
    Ki_PP1_I1       = 1.0e-3          # (uM) Ki for PP1 inhibition by I-1
    kf_PP1_I1       = 1               # (uM) Ki for PP1 inhibition by I-1
    PP1tot                      = PP1_PLBtot      # PP1tot = 0.89  # (uM) total phosphatase 1
    kr_PP1_I1                   = Ki_PP1_I1       # (uM) Ki for PP1 inhibition by I-1
    
    I1 = I1tot - I1ptot
    PP1 = PP1tot - I1p_PP1
    I1p = I1ptot - I1p_PP1
    I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1)
    I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot)

    PP1_eqs = [
        D(I1p_PP1) ~ (kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1),
        D(I1ptot) ~ (I1_phosph - I1_dephosph)
    ]


    ## PLB module
    PLBtot = PLBtotBA   # [uM]
    k_PKA_PLB = 54e-3   # [1/ms]
    Km_PKA_PLB = 21     # [uM]
    k_PP1_PLB = 8.5e-3  # [1/ms]
    Km_PP1_PLB = 7.0    # [uM]
    
    PLB = PLBtot - PLBp
    PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB)
    PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp)
    
    PLB_eqs = [
        D(PLBp) ~ (PLB_phosph - PLB_dephosph)
    ]


    ## PLM module (included 09/18/12) MOUSE
    PLMtot = PLMtotBA   # [uM]
    k_PKA_PLM = 54e-3   # [1/ms]
    Km_PKA_PLM = 21     # [uM]
    k_PP1_PLM = 8.5e-3  # [1/ms]
    Km_PP1_PLM = 7.0    # [uM]

    PLM = PLMtot - PLMp
    PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM)
    PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp)
    
    PLM_eqs = [
        D(PLMp) ~ (PLM_phosph - PLM_dephosph)
    ]

    ## LCC module
    PKAIItot = 0.059        # (uM) total type 2 PKA # MOUSE
    LCCtot = LCCtotBA       # [uM]
    PKACII_LCCtot = 0.025   # [uM]
    PP1_LCC = 0.025         # [uM]
    PP2A_LCC = 0.025        # [uM]
    k_PKA_LCC = 54e-3       # [1/ms]
    Km_PKA_LCC = 21         # [uM]
    k_PP1_LCC = 8.52e-3     # [1/ms] RABBIT, MOUSE
    Km_PP1_LCC = 3          # [uM]
    k_PP2A_LCC = 10.1e-3    # [1/ms]
    Km_PP2A_LCC = 3         # [uM]

    PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII
    LCCa = LCCtot - LCCap
    LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa)
    LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap)
    LCCb = LCCtot - LCCbp
    LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb)
    LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp)
    
    LCC_eqs = [
        D(LCCap) ~ (LCCa_phosph - LCCa_dephosph),
        D(LCCbp) ~ (LCCb_phosph - LCCb_dephosph)
    ]

    ## RyR module (not included in Yang-Saucerman)
    PKAIIryrtot = 0.034         # [uM]
    PP1ryr = 0.034              # [uM]
    PP2Aryr = 0.034             # [uM]
    kcat_pka_ryr = 54e-3        # [1/ms]
    Km_pka_ryr = 21             # [uM]
    kcat_pp1_ryr = 8.52e-3      # [1/ms]
    Km_pp1_ryr = 7              # [uM]
    kcat_pp2a_ryr = 10.1e-3     # [1/ms]
    Km_pp2a_ryr = 4.1           # [uM]
    RyRtot = 0.135              # [uM]
    
    PKACryr = (PKAIIryrtot/PKAIItot)*PKACII
    RyR = RyRtot-RyRp
    RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR)
    RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp)
    RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp)

    RyR_eqs = [
        D(RyRp) ~ (RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A)
    ]


    ## TnI module
    PP2A_TnI = 0.67         # [uM]
    k_PKA_TnI = 54e-3       # [1/ms]
    Km_PKA_TnI = 21         # [uM]
    k_PP2A_TnI = 10.1e-3    # [1/ms]
    Km_PP2A_TnI = 4.1       # [uM]
    TnItot = 70             # [uM] 

    TnI = TnItot - TnIp
    TnI_phosph = k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI)
    TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp)
    
    TnI_eqs = [
        D(TnIp) ~ (TnI_phosph - TnI_dephosph)
    ]
    

    ## Iks module (not present in mouse)

    IKstot = 0.025
    # p(79) = IKstot  # Iks_tot       [uM]
    # p(80) = 0.025   # Yotiao_tot    [uM]
    # p(81) = 0.1e-3  # K_yotiao      [uM] ** apply G589D mutation here **
    # p(82) = 0.025   # PKAII_ikstot  [uM]
    # p(83) = 0.025   # PP1_ikstot    [uM]
    # p(84) = 54      # k_pka_iks     [1/sec]
    # p(85) = 21      # Km_pka_iks    [uM]
    # p(86) = 8.52    # k_pp1_iks     [1/sec]
    # p(87) = 7       # Km_pp1_iks    [uM]
    # 
    # IksYot = y(27)*y(28)/p(81)            # [uM]
    # ydot(27) = p(79) - IksYot - y(27)     # [uM]
    # ydot(28) = p(80) - IksYot - y(28)     # [uM]
    # PKACiks = (IksYot/p(79))*(p(82)/p(34))*y(18) 
    # PP1iks = (IksYot/p(79))*p(83) 
    # Iks = p(79)-y(29) 
    # IKS_PHOSPH = p(40)*p(84)*PKACiks*Iks/(p(85)+p(40)*Iks) 
    # IKS_DEPHOSPH = p(40)*p(86)*PP1iks*y(29)/(p(87)+p(40)*y(29)) 
    # ydot(29) = IKS_PHOSPH-IKS_DEPHOSPH 
    # end Iks module

    IKs_eqs = [
        D(KS79) ~ 0,  # ydot(27) not ODE
        D(KS80) ~ 0,  # ydot(28) not ODE
        D(KSp) ~ 0  # ydot(29)
    ]
    

    ## CFTR module (included 04/30/10)

    ICFTRtot = 0.025  # p(88) = ICFTRtot   # ICFTR_tot      [uM]
    # PKAII_CFTRtot = 0.025   #p(89) = 0.025   # PKAII_CFTRtot [uM]
    # PP1_CFTRtot = 0.025   #p(90) = 0.025   # PP1_CFTRtot   [uM]
    # k_pka_CFTR = 54e-3      #p(91) = 54      # k_pka_CFTR    [1/ms]
    # Km_pka_CFTR = 8.5     #p(92) = 8.5     # Km_pka_CFTR   [uM]
    # k_pp1_CFTR = 8.52e-3    #p(93) = 8.52    # k_pp1_CFTR    [1/ms]
    # Km_pp1_CFTR = 7       #p(94) = 7       # Km_pp1_CFTR   [uM]
    # 
    # CFTRn = ICFTRtot - CFTRp   # Non-phos = tot - phos
    # PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*PKACII     # (PKACFTRtot/PKAIItot)*PKAIIact
    # CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn) 
    # CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*CFTRp/(Km_pp1_CFTR + epsilon*CFTRp)
    
    CFTR_eqs = [
        D(CFTRp) ~ 0  #CFTRphos - CFTRdephos  # ydot(30)
    ]
    


    ## Ikur module (included 04/10/12) MOUSE
    PKAII_KURtot = 0.025    # [uM]
    PP1_KURtot = 0.025      # [uM]
    k_pka_KUR = 54e-3       # [1/ms]
    Km_pka_KUR = 21         # [uM]
    k_pp1_KUR = 8.52e-3     # [1/ms]
    Km_pp1_KUR = 7          # [uM]
    IKurtot = 0.025         # [uM]

    KURn = IKurtot - KURp   # Non-phos = tot - phos
    PKAC_KUR = (PKAII_KURtot/PKAIItot)*PKACII     # (PKA_KURtot/PKAIItot)*PKAIIact
    KURphos = epsilon*KURn*PKAC_KUR*k_pka_KUR/(Km_pka_KUR+epsilon*KURn) 
    KURdephos = PP1_KURtot*k_pp1_KUR*epsilon*KURp/(Km_pp1_KUR+epsilon*KURp) 
    
    Ikur_eqs = [
        D(KURp) ~ (KURphos - KURdephos)
    ]

    #eqs = append!(bar_eqs, cAMP_eqs, PKA_eqs, PP1_eqs, PLB_eqs, PLM_eqs, LCC_eqs, RyR_eqs, TnI_eqs, IKs_eqs, CFTR_eqs, Ikur_eqs)
    #@named sys = ODESystem(eqs, t)
    #return sys  
    return vcat(bar_eqs, cAMP_eqs, PKA_eqs, PP1_eqs, PLB_eqs, PLM_eqs, LCC_eqs, RyR_eqs, TnI_eqs, IKs_eqs, CFTR_eqs, Ikur_eqs)
end
