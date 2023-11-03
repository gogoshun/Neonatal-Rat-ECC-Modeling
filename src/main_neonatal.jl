#=
    This file integrate all ODE files.
=#

include("camsl_ODEfile(M).jl")
include("camcyt_ODEfile(M).jl")
include("camkii_ODEfile(M).jl")
include("ecc_neonatal.jl")
include("bar_new.jl")

using ModelingToolkit
using DifferentialEquations
using Plots
using Sundials


#eq_camdyad = get_camdyad_equations()
eq_camsl = get_camsl_equations()
eq_camcyt = get_camcyt_equations()
eq_camkii = get_camkii_equations()
eq_ecc =  get_ecc_equations()
eq_bar = get_bar_equations()


@named sys = ODESystem(vcat(eq_camdyad, eq_camcyt, eq_camsl, eq_camkii, eq_ecc, eq_bar))


sys = structural_simplify(sys)

#print(states(sys))

@unpack Cai, i_CaJSR, CaNSR, V, Nai, Ki, i_b, i_g, i_d, i_f, i_fca, i_y, i_r, i_s, i_sslow,
        i_Nam, i_Nah, i_Naj, i_PO1, i_PO2, i_PC2, i_nKs, i_CK1, i_CK2, i_OK, i_IK, Cai_sub_SR, Cai_sub_SL,                          # ecc_ODEfile
        CaM_sl, Ca2CaM_sl, Ca4CaM_sl, CaMB_sl, Ca2CaMB_sl, Ca4CaMB_sl, Pb2_sl, Pb_sl, 
        Pt_sl, Pt2_sl, Pa_sl, Ca4CaN_sl, CaMCa4CaN_sl, Ca2CaMCa4CaN_sl, Ca4CaMCa4CaN_sl,                                            # camsl_ODEfile
        CaM_cyt, Ca2CaM_cyt, Ca4CaM_cyt, CaMB_cyt, Ca2CaMB_cyt, Ca4CaMB_cyt, Pb2_cyt, Pb_cyt, 
        Pt_cyt, Pt2_cyt, Pa_cyt, Ca4CaN_cyt, CaMCa4CaN_cyt, Ca2CaMCa4CaN_cyt, Ca4CaMCa4CaN_cyt,                                     # camcyt_ODEfile
        JCaDyad, JCaCyt, JCaSL,
        LCC_PKAp, LCC_CKdyadp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp,                                                              # camkii_ODEfile
        LR, LRG, RG, b1AR_S464, b1AR_S301, GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, PDEp, cAMPtot, RC_I, RCcAMP_I, 
        RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII,                         # bar_ODEfile 
        PKACII_PKI, I1p_PP1, I1ptot, PLBp, PLMp, LCCap, LCCbp, RyRp, TnIp, KURp, KS79, KS80, KSp, CFTRp = sys
        

oprob = ODEProblem(sys, 
        [Cai[1] => 0.22245, Cai[2] => 0.22275, Cai[3] => 0.22302, Cai[4] => 0.22328, Cai[5] => 0.22352, Cai[6] => 0.22375, 
        Cai[7] => 0.22395, Cai[8] => 0.22413, Cai[9] => 0.2243, Cai[10] => 0.22446, Cai[11] => 0.22459, Cai[12] => 0.22471, 
        Cai[13] => 0.22482, Cai[14] => 0.22491, Cai[15] => 0.22499, Cai[16] => 0.22505, Cai[17] => 0.2251, Cai[18] => 0.22514, 
        Cai[19] => 0.22517, Cai[20] => 0.22518, Cai[21] => 0.22518, Cai[22] => 0.22517, Cai[23] => 0.22515, Cai[24] => 0.22511, 
        Cai[25] => 0.22507, Cai[26] => 0.22501, Cai[27] => 0.22495, Cai[28] => 0.22487, Cai[29] => 0.22478, Cai[30] => 0.22469,
        Cai[31] => 0.22458, Cai[32] => 0.22447, Cai[33] => 0.22447, Cai[34] => 0.22421, Cai[35] => 0.22407, Cai[36] => 0.22392, 
        Cai[37] => 0.22376, Cai[38] => 0.22359, Cai[39] => 0.22341, Cai[40] => 0.22323, Cai[41] => 0.22304, Cai[42] => 0.22284, 
        Cai[43] => 0.22263, Cai[44] => 0.22241, Cai[45] => 0.22219,
        i_CaJSR => 778.1041, CaNSR => CaNSR, V => -69.44916, Nai => -69.44916, Ki => 150958.8, i_b => 0.0027, i_g => 0.63674, 
        i_d => 0.63674, i_f => 0.99862, i_fca => 0.99578, i_y => 0.99578, i_r => 0.00672, i_s => 0.97534, i_sslow => 0.22548,
        i_Nam => 0.02272, i_Nah => 0.24062, i_Naj => 0.20256, i_PO1 => 0.20256, i_PO2 => 0, i_PC2 => 0, i_nKs => 0.09345, 
        i_CK1 => 0.00174, i_CK2 => 0.01054, i_OK => 0.27885, i_IK =>0.08076, Cai_sub_SR => 0.22245, Cai_sub_SL => 0.22219, 
        CaM_sl => 4.42e-2, Ca2CaM_sl => 7.34e-5, Ca4CaM_sl => 8.89e-9, 
        CaMB_sl => 2.44, Ca2CaMB_sl => 11.86, Ca4CaMB_sl => 4.38e-4, Pb2_sl => 1.47e-5, Pb_sl => 6.31e-6, Pt_sl => 6.60e-8, 
        Pt2_sl => 7.37e-13, Pa_sl => 4.37e-9, Ca4CaN_sl => 5.22e-4, CaMCa4CaN_sl => 1.98e-6, Ca2CaMCa4CaN_sl => 5.02e-6, 
        Ca4CaMCa4CaN_sl => 1.43e-3, CaM_cyt => 4.4e-2, Ca2CaM_cyt => 4.11e-5, Ca4CaM_cyt => 6.17e-10, CaMB_cyt => 4.179, 
        Ca2CaMB_cyt => 1.11, Ca4CaMB_cyt => 1.61e-5, Pb2_cyt => 8.23e-6, Pb_cyt => 4.15e-8, Pt_cyt => 2.29e-13, 
        Pt2_cyt => 2.47e-18, Pa_cyt => 1.53e-14, Ca4CaN_cyt => 1.17e-4, CaMCa4CaN_cyt => 4.39e-7, Ca2CaMCa4CaN_cyt => 1.59e-7, 
        Ca4CaMCa4CaN_cyt => 1.49e-6, LCC_PKAp => 16.454, LCC_CKdyadp =>16.934 , RyR2809p => 297.36, RyR2815p => 76.985, 
        PLBT17p => 0.614, LCC_CKslp => 8.66e-6, LR => -6.7e-36, LRG => 2.46e-34, RG => 4.8e-4, b1AR_S464 => 5.97e-35, 
        b1AR_S301 => 6.48e-4, GsaGTPtot => 9.6e-3, GsaGDP => 6.21e-4, Gsby => 0.01, AC_GsaGTP => 1.42e-3, PDEp => 2.22e-3, 
        cAMPtot => 1.023, RC_I => 0.804, RCcAMP_I => 0.142, RCcAMPcAMP_I => 4.48e-3, RcAMPcAMP_I => 0.229, PKACI => 8.55e-2, 
        PKACI_PKI => 0.144, RC_II => 0.051, RCcAMP_II => 8.99e-3, RCcAMPcAMP_II => 2.84e-4, RcAMPcAMP_II => 5.77e-2, 
        PKACII => 2.15e-2, PKACII_PKI => 3.62e-2, I1p_PP1 => 7.27e-2, I1ptot => 7.28e-2, PLBp => 8.454, PLMp => 5.6, 
        LCCap => 5.49e-3, LCCbp => 6.27e-3, RyRp => 2.76e-2, TnIp => 4.389, KURp => 1.09e-2, KS79 => 1.53e-3, KS80 => 1.53e-3, 
        KSp => 1.84e-3, CFTRp => 4.06e-3,JCaDyad => 0.0, JCaCyt => 0.0, JCaSL => 0.0], 
        (0.0, 100000.0))

sol = solve(oprob, Rodas5(), tstops = 0:1000:100000, abstol = 1e-10, reltol = 1e-10)
#, abstol = 1e-9, reltol = 1e-9,plotdensity=100
# Solver translated from MATLAB: 
# ode15s/vode –> QNDF() or FBDF(), though in many cases Rodas4(), KenCarp4(), TRBDF2(), or RadauIIA5() are more efficient
length(sol.t)

plot(sol, idxs=Cai, linewidth=1.5, title="Calcium Transient",fmt=:png, xlabel="Time(ms)", ylabel="[Ca]i(mM)",ylim=(0.00005,0.0007),xlim=(0,10.15e4),label="ISO=0.0")

plot(sol, idxs=Vm, linewidth=1.5, title="Action Potential", xlabel="Time(ms)", ylabel="Voltage (mV)",ylim=(-90,60),xlim=(0,10.15e4),label="ISO=0.0")

plot(sol, idxs=Vm, linewidth=3, title="Vm Transient", xlabel="Time(ms)", ylabel="[Na]i(mM)")

plot(sol, idxs=[LCC_CKslp], linewidth=3, xlabel="Time(ms)", ylabel="JCaDyad",xlim=(0,10.15e4))

plot(sol, idxs=[Pb_dyad+Pb2_dyad+Pt_dyad+Pt2_dyad+Pa_dyad]/100, linewidth=1.5, title="CaMKII Activity", xlabel="Time(ms)", ylabel="dCaMKII(mM)",label="with ISO",ylim=(8,10),xlim=(0,5.15e4))
