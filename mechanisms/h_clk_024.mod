TITLE I-h channel from Magee 1998 for distal dendrites
: default values are for dendrites and low Na
: plus leakage, M.Migliore Mar 2010

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
        ehd  		(mV)        
	celsius 	(degC)
	ghdbar=.0001 	(mho/cm2)
        vhalfl=-77.457  	(mV)
        vhalft=-70.242  	(mV)
        a0t=0.0046744   	(/ms)
        zetal=3.4798    	(1)
        zetat=7.3336    	(1)
        gmt=0.14472       (1)
	q10=4.5
	qtl=1
	clk=0.24064
	elk = -70 (mV)
	sh=0
}


NEURON {
	THREADSAFE SUFFIX hdpas
	NONSPECIFIC_CURRENT i
	NONSPECIFIC_CURRENT lk
        RANGE ghdbar, elk, glk, sh 
        GLOBAL linf, taul, clk, vhalfl
}


STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
	lk (mA/cm2)
        linf      
        taul
        ghd
	glk
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	i = ghd*(v-ehd)
	lk = clk*ghdbar*(v-elk)
}


FUNCTION alpl(v(mV)) {
  alpl = exp(0.0378*zetal*(v-vhalfl-sh)) 
}

FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft-sh)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft-sh)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1+ alpl(v))
        taul = bett(v)/(qtl*qt*a0t*(1+a))
}




