COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak conductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 is very small compared to tau1, this is an alphasynapse with time constant tau2.
If tau1/tau2 is very small, this is single exponential decay with time constant tau2.

The factor is evaluated in the initial block 
such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS mygaba
	RANGE tau1, tau2, e, i, A, B, aW, decay, decay_add, deltat, deltat_inv, tlast
	RANGE tau_d, tau_a, tau_decay, tau_da
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau_decay = 1e3
	tau_da = 10
	tau_d = 10
	tau_a = 10
	
	tau1 = 0.1 (ms) <1e-9,1e9>
	tau2 = 3 (ms) <1e-9,1e9>
	e=0	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	spike
	deltat
	deltat_inv
	tlast
	rise_add
}

STATE {
	A (uS)
	B (uS)
	aW
	spikeadder
	rise
	decay
	decay_add
	
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	tlast = 0
	spike = 0
	decay_add = 0
	rise_add = 0
	spikeadder = 0
	decay = 0
	rise = 0
	
	
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	decay' = -decay/tau_decay + decay_add*40
	decay_add' = -decay_add/tau_da
	
	aW = 1/(decay+1)
	
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	deltat = t - tlast
	deltat_inv = 1/deltat

	decay_add = deltat_inv*0.1

	
	if (deltat > 100) {
	  deltat = 100
	}
	
	tlast = t
  
	A = A + aW*weight*factor
	B = B + aW*weight*factor
}

