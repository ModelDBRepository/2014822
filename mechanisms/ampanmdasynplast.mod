NEURON {
    POINT_PROCESS AMPANMDASynPlast
    USEION na WRITE ina
    USEION ca WRITE ica

    : Parameters
    RANGE gAMPAbar, e, AMPA_tau1, AMPA_tau2, gNMDAbar, gNR2Abar, gNR2Bbar
    RANGE eNa, eCa
    RANGE mg
    RANGE Alpha_nr2a, Alpha_nr2b, Beta_nr2a, Beta_nr2b
    RANGE v_thresh1, v_thresh2, v_trace1_thresh2, v_trace1_tau, v_trace2_tau
    RANGE A_ltp, A_ltd
    RANGE wmax, wmin
    RANGE hill_coef_ltp, hill_coef_ltd, hill_midpoint_ltp, hill_midpoint_ltp2, hill_midpoint_ltd

    : AMPA
    RANGE g_ampa, i_ampa, Cmax, Alpha_ampa, Beta_ampa

    : NMDA
    RANGE g_nr2a, g_nr2b, g_nmda, i_nmda, g_nr2a_base, g_nr2b_base, g_nmda_base
    RANGE g_nmda_LTP, g_nmda_LTD
    
    
    : SynPlast
    RANGE g_nmda_trace_ltp_tau, g_nmda_trace_ltd_tau, X_trace_tau
    RANGE g_nmda_trace_ltp, g_nmda_trace_ltd, X_trace, X, X_max
    RANGE v_threshed1, v_threshed2
    RANGE ltp_part, ltd_part, ltp_mult, ltd_mult, w_change, ltp_nr2a_part, ltp_nr2b_part
    RANGE v_trace1_threshed1, v_trace2_threshed2
    RANGE w_ltp, w_ltd, weight, hilleq_ltp, hilleq_ltp2, hilleq_ltd
    RANGE hill_nmda_ltd, hill_nmda_ltp
    RANGE moving_threshold_hill_ltd_multiplier, moving_threshold_hill_ltd_tau, moving_threshold_hill_ltd
    RANGE moving_threshold_hill_ltp_multiplier, moving_threshold_hill_ltp_tau, moving_threshold_hill_ltp
    RANGE moving_threshold_hill_ltp2_multiplier, moving_threshold_hill_ltp2_tau, moving_threshold_hill_ltp2
    RANGE nmda_thresh, nmda_threshed
    RANGE g_nmda_LTP, g_nmda_LTD
    : Misc

    : Debug


    GLOBAL Cdur, Rinf_nr2a, Rinf_nr2b, Rtau_nr2a, Rtau_nr2b
    NONSPECIFIC_CURRENT i
}

UNITS {
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (nS) = (nanosiemens)
    (umho) = (micromho)
    (mM) = (milli/liter)
}

PARAMETER {
    : ----------------------------------- AMPA
    gAMPAbar = 1    (uS)

    : - Current
    e = 0	(mV)

    : - Conductance
    AMPA_tau1 = .1 (ms) <1e-9,1e9>
    AMPA_tau2 = 10 (ms) <1e-9,1e9>

    : ----------------------------------- Destexhe AMPA
    Alpha_ampa = 1.1		(/ms)	: forward (binding) rate
    Beta_ampa = 0.19		(/ms)	: backward (unbinding) rate
    Cmax = 1

    : ----------------------------------- NMDA
    gNMDAbar = 1    (uS)
    gNR2Abar = 1    (uS)
    gNR2Bbar = 1    (uS)

    : - Current
    e_nmda = 0	(mV)

    : - Conductance
    mg = 1    (mM)		: external magnesium concentration
    Cdur = 1		(ms)	: transmitter duration (rising phase)
    Alpha_nr2a = 0.35		(/ms)	: forward (binding) rate
    Alpha_nr2b = 0.35		(/ms)	: forward (binding) rate
    Beta_nr2a = 0.015		(/ms)	: backward (unbinding) rate
    Beta_nr2b = 0.015		(/ms)	: backward (unbinding) rate


    : ----------------------------------- SynPlast
    v_thresh1 = -25.0   (mV)
    v_thresh2 = -58.0   (mV)
    v_trace1_thresh2 = 0.1
    v_trace1_tau = 200  (ms)
    v_trace2_tau = 40   (ms)
    g_nmda_trace_ltp_tau = 50   (ms)
    g_nmda_trace_ltd_tau = 50   (ms)
    X_trace_tau = 2 (ms)
    X_max = 1

    A_ltp = 5.0e-2
    A_ltd = 1.0e-4

    wmax = 2
    wmin = 0
    nmda_thresh = 0.01      (uS)

    : ----------------------------------- Misc

    hill_coef_ltp = 4
    hill_midpoint_ltp = 0.2
    hill_midpoint_ltp2 = 0.2

    hill_coef_ltd = 2
    hill_midpoint_ltd = 0.8

    moving_threshold_hill_ltd_multiplier = 1.0
    moving_threshold_hill_ltd_tau = 10.0    (ms)

    moving_threshold_hill_ltp_multiplier = 1.0
    moving_threshold_hill_ltp_tau = 10.0    (ms)

    moving_threshold_hill_ltp2_multiplier = 1.0
    moving_threshold_hill_ltp2_tau = 10.0   (ms)

    : ----------------------------------- Debug
    
}

ASSIGNED {
    : AMPA
    i_ampa (nA)
    : factor

    g_ampa  (uS)
    : Destexhe AMPA

    
    synon_ampa
    Rinf_ampa
    Rtau_ampa       (ms)

    : NMDA
    i_nmda (nA)
    Rinf_nr2a				: steady state channels open
    Rinf_nr2b				: steady state channels open
    Rtau_nr2a		(ms)		: time constant of channel binding
    Rtau_nr2b		(ms)		: time constant of channel binding

    synon
    mg_mult

    g_nr2a      (uS)
    g_nr2b      (uS)
    g_nr2a_base_wbar
    g_nr2b_base_wbar
    g_nr2a_base
    g_nr2b_base
    g_nmda      (uS)
    g_nmda_base (uS)
    g_nmda_LTP  (uS)
    g_nmda_LTD  (uS)
    
    ina     (nA) : (mA/cm2)  : Na current density
    ica     (nA) : (mA/cm2)  : Ca current density

    : SynPlast
    v_threshed1  (mV)
    v_threshed2  (mV)


    hilleq_ltp
    hilleq_ltp2
    hilleq_ltd
    hill_nmda_ltp   (uS)
    hill_nmda_ltd   (uS)

    ltp_nr2a_part
    ltp_nr2b_part
    ltp_part
    ltd_part

    ltp_mult
    ltd_mult

    w_change
    X

    : Misc
    v (mV)
    i (nA)


    : Debug
}

STATE {
    : AMPA
    A (uS)
    B (uS)

    : Destexhe AMPA
    Ron_ampa
    Roff_ampa

    : NMDA
    Ron_nr2a
    Ron_nr2b
    Roff_nr2a
    Roff_nr2b

    : SynPlast
    v_trace1_threshed1      (mV)
    v_trace2_threshed2      (mV)
    g_nmda_trace_ltp        (uS)
    g_nmda_trace_ltd        (uS)
    X_trace

    moving_threshold_hill_ltd
    moving_threshold_hill_ltp
    moving_threshold_hill_ltp2

    nmda_threshed       (uS)

    w_ltp
    w_ltd

    weight
    : Misc

    : Debug
}

INITIAL {
    LOCAL tp
    : AMPA
    Rinf_ampa = Alpha_ampa / (Alpha_ampa + Beta_ampa)
    Rtau_ampa = 1 / (Alpha_ampa + Beta_ampa)

    Ron_ampa = 0
    Roff_ampa = 0
    synon_ampa = 0

    : NMDA
    Rinf_nr2a = Alpha_nr2a / (Alpha_nr2a + Beta_nr2a)
    Rinf_nr2b = Alpha_nr2b / (Alpha_nr2b + Beta_nr2b)
    Rtau_nr2a = 1 / (Alpha_nr2a + Beta_nr2a)
    Rtau_nr2b = 1 / (Alpha_nr2b + Beta_nr2b)

    Ron_nr2a = 0
    Ron_nr2b = 0
    Roff_nr2a = 0
    Roff_nr2b = 0
    
    synon = 0

    : SynPlast
    v_trace1_threshed1 = 0
    v_trace2_threshed2 = 0

    nmda_threshed = 0
    w_ltp = 0
    w_ltd = 0

    v_threshed1 = 0
    v_threshed2 = 0
    g_nmda_trace_ltp = 0
    g_nmda_trace_ltd = 0
    X_trace = 0

    moving_threshold_hill_ltd = 0
    moving_threshold_hill_ltp = 0
    moving_threshold_hill_ltp2 = 0

    hilleq_ltp = 0
    hilleq_ltp2 = 0
    hilleq_ltd = 0
    hill_nmda_ltp = 0
    hill_nmda_ltd = 0

    ltp_nr2a_part = 0
    ltp_nr2b_part = 0
    ltp_part = 0
    ltd_part = 0

    ltp_mult = 0
    ltd_mult = 0

    w_change = 0
    weight = 1
    X = 0

    : Misc

    : Debug

}

BREAKPOINT {
    SOLVE state METHOD cnexp

    if (weight < 0) {
      weight = 0
    }
    : Conductance
    g_ampa = (Ron_ampa + Roff_ampa) * gAMPAbar * weight

    : - NMDA
    mg_mult = mgblock(v)
    g_nr2a_base = mg_mult*(Ron_nr2a + Roff_nr2a)
    g_nr2b_base = mg_mult*(Ron_nr2b + Roff_nr2b)

    g_nr2a_base_wbar = g_nr2a_base * gNR2Abar
    g_nr2b_base_wbar = g_nr2b_base * gNR2Bbar
    g_nr2a = g_nr2a_base * gNMDAbar
    g_nr2b = g_nr2b_base * gNMDAbar
    g_nmda_base = (g_nr2a_base_wbar + g_nr2b_base_wbar)
    g_nmda = g_nmda_base * gNMDAbar
   

    g_nmda_LTP = (0.2 * g_nr2a_base_wbar + 0.8 * g_nr2b_base_wbar)
    g_nmda_LTD = (0.8 * g_nr2a_base_wbar + 0.2 * g_nr2b_base_wbar)





    : Current
    : - NMDA
    ina = 0.94 * g_nmda * (v - e_nmda)   : current in mA/cm2 
    ica = 0.06 * g_nmda * (v - e_nmda)   : current in mA/cm2

    if (v>=e_nmda) {
       ica = 0
    }

    i_ampa = g_ampa * (v - e)
    i_nmda = (ina + ica)
    i = i_ampa + i_nmda
}

DERIVATIVE state {
    : AMPA
    Ron_ampa' = (synon_ampa * Rinf_ampa - Ron_ampa) / Rtau_ampa
    Roff_ampa' = -Beta_ampa * Roff_ampa

    : NMDA
    Ron_nr2a' = (synon*Rinf_nr2a - Ron_nr2a)/Rtau_nr2a
    Ron_nr2b' = (synon*Rinf_nr2b - Ron_nr2b)/Rtau_nr2b
    Roff_nr2a' = -Beta_nr2a*Roff_nr2a
    Roff_nr2b' = -Beta_nr2b*Roff_nr2b

    : Plasticity
    : - Vmem LTP
    v_threshed1 = v - v_thresh1
    if (v_threshed1 < 0) {
      v_threshed1 = 0
    }

    v_trace1_threshed1' = (v_threshed1 - v_trace1_threshed1) / v_trace1_tau

    : - Vmem LTD
    v_threshed2 = v - v_thresh2
    if (v_threshed2 < 0) {
      v_threshed2 = 0
    }

    v_trace2_threshed2' = (v_threshed2 - v_trace2_threshed2) / v_trace2_tau
    
    : - NMDA
    nmda_threshed = g_nmda_base - nmda_thresh
    if (nmda_threshed < 0) {
      nmda_threshed = 0
    }

    g_nmda_trace_ltp' = (g_nmda_LTP - g_nmda_trace_ltp) / g_nmda_trace_ltp_tau
    g_nmda_trace_ltd' = (g_nmda_LTD - g_nmda_trace_ltd) / g_nmda_trace_ltd_tau
    
    : - Dirac
    X_trace' = (X - X_trace) / X_trace_tau
    

    : - Hill LTP
    hill_nmda_ltp = (g_nmda_trace_ltp ^ hill_coef_ltp)
    hilleq_ltp = hill_nmda_ltp / ((hill_midpoint_ltp) ^ hill_coef_ltp + hill_nmda_ltp) - hill_nmda_ltp / ((hill_midpoint_ltp2) ^ hill_coef_ltp + hill_nmda_ltp) - moving_threshold_hill_ltp
    if (hilleq_ltp < 0) {
        hilleq_ltp = 0
    }

    hilleq_ltp2 = hill_nmda_ltp / ((hill_midpoint_ltp2) ^ hill_coef_ltp + hill_nmda_ltp) - moving_threshold_hill_ltp
    if (hilleq_ltp2 < 0) {
        hilleq_ltp2 = 0
    }

    : - Hill LTD
    hill_nmda_ltd = (g_nmda_trace_ltd ^ hill_coef_ltd)
    hilleq_ltd = hill_nmda_ltd / ((hill_midpoint_ltd) ^ hill_coef_ltd + hill_nmda_ltd) - moving_threshold_hill_ltd
    if (hilleq_ltd < 0) {
        hilleq_ltd = 0
    }

    : - Moving Threshold Hill
    moving_threshold_hill_ltp' = (hilleq_ltd * moving_threshold_hill_ltp_multiplier - moving_threshold_hill_ltp) / moving_threshold_hill_ltp_tau
    moving_threshold_hill_ltp2' = (hilleq_ltp2 * moving_threshold_hill_ltp2_multiplier - moving_threshold_hill_ltp2) / moving_threshold_hill_ltp2_tau
    moving_threshold_hill_ltd' = (hilleq_ltp * moving_threshold_hill_ltd_multiplier - moving_threshold_hill_ltd) / moving_threshold_hill_ltd_tau
    
    : - Both
    ltp_part = A_ltp * hilleq_ltp * v_trace1_threshed1
    ltd_part = A_ltd * hilleq_ltd * X_trace * v_trace2_threshed2

    ltp_mult = (wmax - weight)
    ltd_mult = (weight - wmin)
    w_change = (ltp_part * ltp_mult - ltd_part * ltd_mult)
    
    w_ltp = ltp_part * ltp_mult
    w_ltd = ltd_part * ltd_mult
    weight' = w_change
    
}

NET_RECEIVE(w, on, nspike, r0_nr2a, r0_nr2b, r0_ampa, t0 (ms)) {
    LOCAL weight_ampa
    weight_ampa = w * 0.1
    : NMDA:
    : flag is an implicit argument of NET_RECEIVE and  normally 0
    if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
        
	    nspike = nspike + 1
        X = X_max
	    if (!on) {
            : AMPA
            r0_ampa = r0_ampa*exp(-Beta_ampa*(t - t0))
		    Ron_ampa = Ron_ampa + r0_ampa
		    Roff_ampa = Roff_ampa - r0_ampa
		    synon_ampa = synon_ampa + weight_ampa


            : NMDA

		    r0_nr2a = r0_nr2a*exp(-Beta_nr2a*(t - t0))
		    r0_nr2b = r0_nr2b*exp(-Beta_nr2b*(t - t0))
		    
		    on = 1
		    t0 = t
		    synon = synon + w
		    
		    Ron_nr2a = Ron_nr2a + r0_nr2a
		    Roff_nr2a = Roff_nr2a - r0_nr2a
		    
		    Ron_nr2b = Ron_nr2b + r0_nr2b
		    Roff_nr2b = Roff_nr2b - r0_nr2b
	    }
	    : come again in Cdur with flag = current value of nspike
	    net_send(Cdur, nspike)
    }
    if (flag == nspike) { : if this associated with last spike then turn off
        
        : AMPA
	    r0_ampa = weight_ampa*Rinf_ampa + (r0_ampa - weight_ampa*Rinf_ampa)*exp(-(t - t0)/Rtau_ampa)
	    Ron_ampa = Ron_ampa - r0_ampa
	    Roff_ampa = Roff_ampa + r0_ampa
	    synon_ampa = synon_ampa - weight_ampa

        : NMDA

	    r0_nr2a = w*Rinf_nr2a + (r0_nr2a - w*Rinf_nr2a)*exp(-(t - t0)/Rtau_nr2a)
	    r0_nr2b = w*Rinf_nr2b + (r0_nr2b - w*Rinf_nr2b)*exp(-(t - t0)/Rtau_nr2b)
	    synon = synon - w
	    Ron_nr2a = Ron_nr2a - r0_nr2a
	    Ron_nr2b = Ron_nr2b - r0_nr2b
	    Roff_nr2a = Roff_nr2a + r0_nr2a
	    Roff_nr2b = Roff_nr2b + r0_nr2b

	    t0 = t
	    on = 0
        X = 0

    }
}

FUNCTION mgblock(v(mV)) {
    TABLE 
    DEPEND mg
    FROM -140 TO 80 WITH 1000

    : from Jahr & Stevens
    mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

