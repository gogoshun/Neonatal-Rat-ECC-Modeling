#=
    This function computes the CaMKII-dependent phosphorylation profiles for
    LTCCs (dyadic and subsarcolemmal), RyRs, and PLB.
=#
using DifferentialEquations
using ModelingToolkit

export get_camkii_equations


function get_camkii_equations()
    @variables t
    D = Differential(t)

    @variables LCC_PKAp(t) LCC_CKdyadp(t) RyR2809p(t) RyR2815p(t) PLBT17p(t) LCC_CKslp(t) 
    
    ## RATE CONSTANTS and KM VALUES

    # L-Type Ca Channel (LTCC) parameters
    k_ckLCC = 0.4                   # [s^-1]
    k_pp1LCC = 0.1103               # [s^-1] 
    k_pkaLCC = 13.5                 # [s^-1] 
    k_pp2aLCC = 10.1                # [s^-1] 

    KmCK_LCC = 12                   # [uM] 
    KmPKA_LCC = 21                  # [uM] 
    KmPP2A_LCC = 47                 # [uM] 
    KmPP1_LCC = 9                   # [uM] 

    # Ryanodine Receptor (RyR) parameters
    k_ckRyR = 0.4                   # [s^-1] 
    k_pkaRyR = 1.35                 # [s^-1] 
    k_pp1RyR = 1.07                 # [s^-1] 
    k_pp2aRyR = 0.481               # [s^-1] 

    # Basal RyR phosphorylation (numbers based on param estimation)
    kb_2809 = 0.51                  # [uM/s] - PKA site
    kb_2815 = 0.35                  # [uM/s] - CaMKII site

    KmCK_RyR = 12                   # [uM] 
    KmPKA_RyR = 21                  # [uM] 
    KmPP1_RyR = 9                   # [uM] 
    KmPP2A_RyR = 47                 # [uM] 

    # Phospholamban (PLB) parameters
    k_ckPLB = 8e-3                  # [s^-1]
    k_pp1PLB = 0.0428               # [s^-1]

    KmCK_PLB = 12
    KmPP1_PLB = 9

    # Okadaic Acid inhibition params (based on Huke/Bers [2008])
    # Want to treat OA as non-competitive inhibitor of PP1 and PP2A
    Ki_OA_PP1 = 0.78                # [uM] - Values from fit
    Ki_OA_PP2A = 0.037              # [uM] - Values from fit

    # Default PKA level
    PKAc = 95.6 * 0.54

    ## Parameters for CaMKII module
    LCCtotDyad = 31.4*0.9      # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    LCCtotSL = 0.0846          # [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
    RyRtot = 382.6             # [uM] - Total RyR (in Dyad)
    PP1_dyad = 95.7            # [uM] - Total dyadic [PP1]
    PP1_SL = 0.57              # [uM] - Total Subsarcolemmal [PP1]
    PP2A_dyad = 95.76          # [uM] - Total dyadic PP2A
    OA = 0                     # [uM] - PP1/PP2A inhibitor Okadaic Acid
    plb_val=106 # MOUSE
    PLBtot = plb_val           # [uM] - Total [PLB] in cytosolic units

    ## OA inhibition term (non-competitive) for PP1 and PP2A
    OA_PP1 = 1/(1 + (OA/Ki_OA_PP1)^3)
    OA_PP2A = 1/(1 + (OA/Ki_OA_PP2A)^3)

    CaMKIItotDyad = 120             # [uM]
    CaMKIItotSL = 120*8.293e-4      # [uM]
    PP1_PLBtot = 0.89               # [uM] - [umol/L cytosol]

    ## ODE EQUATIONS
    # LTCC states (note: PP2A is acting on PKA site and PP1 on CKII site)
    @variables Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t) Pa_dyad(t) I1p_PP1(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t) Pa_sl(t)

    # Variables related to camdyad_ODEfile
    CaMKIIact_Dyad = CaMKIItotDyad .* (Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad)
    CaMKIIact_SL = CaMKIItotSL .* (Pb_sl + Pt_sl + Pt2_sl + Pa_sl)
    PP1_PLB_avail = 1 - I1p_PP1/PP1_PLBtot + 0.081698
    # CaMKII phosphorylation of Dyadic LCCs
    LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp 
    LCCDyad_PHOS = (k_ckLCC*CaMKIIact_Dyad*LCC_CKdyadn)/(KmCK_LCC+LCC_CKdyadn) 
    LCCDyad_DEPHOS = (k_pp1LCC*PP1_dyad*LCC_CKdyadp)/(KmPP1_LCC+LCC_CKdyadp)*OA_PP1 
    LCC_CKsln = LCCtotSL - LCC_CKslp 
    LCCSL_PHOS = (k_ckLCC*CaMKIIact_SL*LCC_CKsln)/(KmCK_LCC+LCC_CKsln) 
    LCCSL_DEPHOS = (k_pp1LCC*PP1_SL*LCC_CKslp)/(KmPP1_LCC+LCC_CKslp)*OA_PP1 
    LCC_PKAn = LCCtotDyad - LCC_PKAp 
    RyR2815n = RyRtot - RyR2815p 
    RyR_BASAL = kb_2815*RyR2815n 
    RyR_PHOS = (k_ckRyR*CaMKIIact_Dyad*RyR2815n)/(KmCK_RyR+RyR2815n) 
    RyR_PP1_DEPHOS = (k_pp1RyR*PP1_dyad*RyR2815p)/(KmPP1_RyR+RyR2815p)*OA_PP1 
    RyR_PP2A_DEPHOS = (k_pp2aRyR*PP2A_dyad*RyR2815p)/(KmPP2A_RyR+RyR2815p)*OA_PP2A 
    RyR2809n = RyRtot - RyR2809p 
    PP1_PLB = PP1_dyad*PP1_PLB_avail  # Inhibitor-1 regulation of PP1_dyad included here
    PLBT17n = PLBtot - PLBT17p 
    PLB_PHOS = (k_ckPLB*PLBT17n*CaMKIIact_Dyad)/(KmCK_PLB+PLBT17n) 
    PLB_DEPHOS = (k_pp1PLB*PP1_PLB*PLBT17p)/(KmPP1_PLB+PLBT17p)*OA_PP1 

    CaMKII_eqs = [
        D(LCC_CKdyadp) ~ (LCCDyad_PHOS - LCCDyad_DEPHOS)*1e-3, # du[2]
        # CaMKII phosphorylation of Sub-sarcolemmal LCCs
        D(LCC_CKslp) ~ (LCCSL_PHOS - LCCSL_DEPHOS)*1e-3, # du[6]
        # PKA phosphorylation (currently unused elsewhere)
        D(LCC_PKAp) ~ ((k_pkaLCC*PKAc*LCC_PKAn)/(KmPKA_LCC+LCC_PKAn) - (k_pp2aLCC*PP2A_dyad*LCC_PKAp)/(KmPP2A_LCC+LCC_PKAp)*OA_PP2A)*1e-3, # du[1]
        # RyR states
        D(RyR2815p) ~ (RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS)*1e-3, # du[4]
        # PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
        D(RyR2809p) ~ (kb_2809*RyR2809n + (k_pkaRyR*PKAc*RyR2809n)/(KmPKA_RyR+RyR2809n) - (k_pp1RyR*PP1_dyad*RyR2809p)/(KmPP1_RyR+RyR2809p)*OA_PP1)*1e-3, # du[3]
        # PLB states
        D(PLBT17p) ~ (PLB_PHOS - PLB_DEPHOS)*1e-3 # du[5]
    ]

    #@named sys = ODESystem(CaMKII_eqs, t)

    return CaMKII_eqs
end
