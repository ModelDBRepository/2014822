COMMENT
Revised 12/15/2000 in light of a personal communication 
from Misha Tsodyks that u is incremented _before_ x is 
converted to y--a point that was not clear in the paper.
If u is incremented _after_ x is converted to y, then 
the first synaptic activation after a long interval of 
silence will produce smaller and smaller postsynaptic 
effect as the length of the silent interval increases, 
eventually becoming vanishingly small.

Implementation of a model of short-term facilitation and depression 
based on the kinetics described in
  Tsodyks et al.
  Synchrony generation in recurrent networks 
  with frequency-dependent synapses
  Journal of Neuroscience 20:RC50:1-5, 2000.
Their mechanism represented synapses as current sources.
The mechanism implemented here uses a conductance change instead.

The basic scheme is

x -------> y    Instantaneous, spike triggered.
                Increment is u*x (see discussion of u below).
                x == fraction of "synaptic resources" that have 
                     "recovered" (fraction of xmtr pool that is 
                     ready for release, or fraction of postsynaptic 
                     channels that are ready to be opened, or some 
                     joint function of these two factors)
                y == fraction of "synaptic resources" that are in the 
                     "active state."  This is proportional to the 
                     number of channels that are open, or the 
                     fraction of max synaptic current that is 
                     being delivered. 
  tau_1
y -------> z    z == fraction of "synaptic resources" that are 
                     in the "inactive state"

  tau_rec
z -------> x

where x + y + z = 1

The active state y is multiplied by a synaptic weight to compute
the actual synaptic conductance (or current, in the original form 
of the model).

In addition, there is a "facilition" term u that 
governs the fraction of x that is converted to y 
on each synaptic activation.

  -------> u    Instantaneous, spike triggered, 
                happens _BEFORE_ x is converted to y.
                Increment is U*(1-u) where U and u both 
                lie in the range 0 - 1.
  tau_facil
u ------->      decay of facilitation

This implementation for NEURON offers the user a parameter 
u0 that has a default value of 0 but can be used to specify 
a nonzero initial value for u.

When tau_facil = 0, u is supposed to equal U.

Note that the synaptic conductance in this mechanism 
has the same kinetics as y, i.e. decays with time 
constant tau_1.

This mechanism can receive multiple streams of 
synaptic input via NetCon objects.  
Each stream keeps track of its own 
weight and activation history.

The printf() statements are for testing purposes only.
ENDCOMMENT


NEURON {
	POINT_PROCESS ltpltd
	RANGE e, i, g2, peso, U, tau_facil, tau_rec, tau_1, Rin, Ase, taum, lambdap, lambdad, gNMDAbar, gAMPAbar
	RANGE u0, Pini, Nini, deltap, deltad, f, nip, nid, gamma, eta, VVini, maxW, gmult, Np, Nd, g_nmda, i_ampa, synweight
	RANGE mg, isynAMPA, ap, ad, mp, md, VV, g_nr2a, g_nr2b, gNR2Abar, gNR2Bbar, Rinf_nr2a, Rtau_nr2a, Rinf_nr2b, Rtau_nr2b, Alpha_nr2a, Beta_nr2a, Alpha_nr2b, Beta_nr2b, g_adj
	GLOBAL Cdur
	NONSPECIFIC_CURRENT i 
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	mg = 1 (mM)
	gNMDAbar = 1
	gAMPAbar = 1
	Cdur	= 1		(ms)	: transmitter duration (rising phase)
	Alpha_nr2a	= 0.35		(/ms)	: forward (binding) rate
	Beta_nr2a	= 0.015		(/ms)	: backward (unbinding) rate
	Alpha_nr2b	= 0.35		(/ms)	: forward (binding) rate
	Beta_nr2b	= 0.015		(/ms)	: backward (unbinding) rate

	: e = -90 mV for inhibitory synapses,
	:     0 mV for excitatory
	e = 0	(mV)
	: tau_1 was the same for inhibitory and excitatory synapses
	: in the models used by T et al.
	tau_1 = 3 (ms) < 1e-9, 1e9 >
	: tau_rec = 100 ms for inhibitory synapses,
	:           800 ms for excitatory
	tau_rec = 50 (ms) < 1e-9, 1e9 >
	: tau_facil = 1000 ms for inhibitory synapses,
	:             0 ms for excitatory
	tau_facil = 200 (ms) < 0, 1e9 >
	: U = 0.05 for inhibitory synapses, 
	:     0.5 for excitatory
	: the (1) is needed for the < 0, 1 > to be effective
	:   in limiting the values of U and u0
	U = 0.36(1) < 0, 1 >
	: initial value for the "facilitation variable"
	u0 = 0 (1) < 0, 1 > 
	:PARAMETERS 
	f = 0.05e-3 
	deltap = 400
	deltad = 400
	gamma = 0.2
	eta = 2e-3
	nip = 0.0987
	nid = 0.07
	lambdap = 1e-3
	lambdad = 2e-3
	mp = 3e-3
	md = 3e-3
	ap = 2
	ad = 0.5
	taum = 40
	Rin=10e7
	Ase=2.5e-7
	Pini=0
	Nini=0
	VVini=0
	g2=43 (umho)
	peso=0.5e-6 
	gmult=1e-3
	}

ASSIGNED {
	g_adj
	synweight
	v (mV)
	i (nA)
	x
	y
	g_nmda
	g_nr2a
	g_nr2b
	i_ampa
	Rinf_nr2a
	Rtau_nr2a (ms)
	Rinf_nr2b
	Rtau_nr2b (ms)
	gNR2Abar
	gNR2Bbar
	synon
	isynAMPA
}

STATE {
	g (umho)
	C
	Np
	Nd
	VV (mV)
	maxW
	Ron_nr2a
	Roff_nr2a
	Ron_nr2b
	Roff_nr2b
}

INITIAL {
	synon = 0
	Rinf_nr2a = Alpha_nr2a / (Alpha_nr2a + Beta_nr2a)
	Rinf_nr2b = Alpha_nr2b / (Alpha_nr2b + Beta_nr2b)
	Rtau_nr2a = 1 / (Alpha_nr2a + Beta_nr2a)
	Rtau_nr2b = 1 / (Alpha_nr2b + Beta_nr2b)
	
	VV=VVini
	g=0
	C =0
	Np=Pini
	Nd=Nini
	maxW = 0
}

	
FUNCTION v2(v) { if (v<-65) { v2 = 0	} else {v2 = v+65}
	}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit

	g_nr2a = 3.4*mgblock(v)*(Ron_nr2a + Roff_nr2a)*1(umho)*gNR2Abar
	g_nr2b = 3.4*mgblock(v)*(Ron_nr2b + Roff_nr2b)*1(umho)*gNR2Bbar
	g_nmda = (g_nr2a + g_nr2b)*gNMDAbar
	
	i_ampa = -g2 * VV * gAMPAbar
	i = i_ampa + g_nmda*(v - e)
	maxW = f*(deltap*Np-deltad*Nd)
	isynAMPA = Ase*g

	synweight = ((Np-Nd)/2 + 1)

	: i = -g2 * VV
}

DERIVATIVE state {

	:printf("%g, %g, %g\n", tau_1, tau_rec, tau_facil)
	:NMDA
	Ron_nr2a' = (synon*Rinf_nr2a - Ron_nr2a)/Rtau_nr2a
	Ron_nr2b' = (synon*Rinf_nr2b - Ron_nr2b)/Rtau_nr2b
	Roff_nr2a' = -Beta_nr2a*Roff_nr2a
	Roff_nr2b' = -Beta_nr2b*Roff_nr2b

	g' = -g/tau_1
	g_adj = g*gmult
	C' = gamma*(VV)-eta*C+peso*(v2(v))
	Np' = nip*C-(lambdap+g_adj*deltap)*Np+((mp*Np*Np)/(ap+Np*Np))
	Nd' = nid*C-(lambdad+g_adj*deltad)*Nd+((md*Nd*Nd)/(ad+Nd*Nd))
	VV' = -(VV/taum)+Rin*Ase*g_adj*((1/taum)+f*(deltap*Np-deltad*Nd))
}



		


:NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) 
NET_RECEIVE(weight (umho), on, nspike, r0_nr2a, r0_nr2b, t0, y, z, u, tsyn (ms)) 
{
INITIAL {
: these are in NET_RECEIVE to be per-stream
	y = 0
	z = 0
:	u = 0
	u = u0
	tsyn = t
: this header will appear once per stream
:printf("t\t t-tsyn\t y\t z\t u\t newu\t g\t dg\t newg\t newy\n")
}
     :flag is an implicit argument of NET_RECEIVE and  normally 0
     if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse

	: first calculate z at event-
	:   based on prior y and z
	z = z*exp(-(t - tsyn)/tau_rec)
	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
	: now calc y at event-
	y = y*exp(-(t - tsyn)/tau_1)

	x = 1-y-z

	: calc u at event--
	if (tau_facil > 0) {
		u = u*exp(-(t - tsyn)/tau_facil)
	} else {
		u = U
	}


	if (tau_facil > 0) {
		state_discontinuity(u, u + U*(1-u))
	}


	state_discontinuity(g, g + weight*x*u)
	state_discontinuity(y, y + x*u)
	
	tsyn = t

    : NMDA:
	    nspike = nspike + 1
	    if (!on) {
		    r0_nr2a = r0_nr2a*exp(-Beta_nr2a*(t - t0))
		    r0_nr2b = r0_nr2b*exp(-Beta_nr2b*(t - t0))
		    
		    on = 1
		    t0 = t
		    synon = synon + weight
		    
		    Ron_nr2a = Ron_nr2a + r0_nr2a
		    Roff_nr2a = Roff_nr2a - r0_nr2a
		    
		    Ron_nr2b = Ron_nr2b + r0_nr2b
		    Roff_nr2b = Roff_nr2b - r0_nr2b
	    }
	    : come again in Cdur with flag = current value of nspike
	    net_send(Cdur, nspike)

    }
    if (flag == nspike) { : if this associated with last spike then turn off
	    r0_nr2a = weight*Rinf_nr2a + (r0_nr2a - weight*Rinf_nr2a)*exp(-(t - t0)/Rtau_nr2a)
	    r0_nr2b = weight*Rinf_nr2b + (r0_nr2b - weight*Rinf_nr2b)*exp(-(t - t0)/Rtau_nr2b)
	    t0 = t
	    synon = synon - weight
	    Ron_nr2a = Ron_nr2a - r0_nr2a
	    Ron_nr2b = Ron_nr2b - r0_nr2b
	    Roff_nr2a = Roff_nr2a + r0_nr2a
	    Roff_nr2b = Roff_nr2b + r0_nr2b
	    on = 0
    }

}



FUNCTION mgblock(v(mV)) {
    TABLE 
    DEPEND mg
    FROM -140 TO 80 WITH 1000

    : from Jahr & Stevens
    mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

