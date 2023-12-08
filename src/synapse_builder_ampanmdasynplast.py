from .synapse_builder import SynapseBuilder
from neuron import h
import os
import pickle
import time
import numpy as np
from pprint import pprint
from pathlib import Path
import pathlib

class AMPANMDASynPlastBuilder(SynapseBuilder):

    epsp_start = 200
    epsp_end = 10000
    after_epsp_wait = 200
    stimulation_start = 1500
    
    root_path = os.getcwd()
    
    def __init__(self):
        super().__init__()
        
        
    def init_object_vars(self):
        super().init_object_vars()
        self.democracy = False
        self.records_synapse = []
        self.records_ec_synapse = []
        self.stimulations = []
        self.clamp_connections = []
        self.network_connections = []
        self.plastic_synapses = []
        self.nonplastic_synapses = []
        self.artificial_cells = []
        self.gaba_synapses = []
        self.weight_factor = 1
        self.plasticity_model_parameters = {}
        self.nr2b_multiplier = 1
        self.gamma_pairings = 3
        self.theta_pairings = 20
        self.gamma_freq = 80
        self.theta_freq = 4
        self.gbarmult = 1.5
        self.nr2b_cure = 1
        self.A_pm_mult = 1
        self.gaba_weight = 0.01
        self.ampamult = 1
        self.nmdamult = 1
        self.nr2amult = 1
        self.nr2bmult = 1
        self.wmax = 2.5
        self.wmin = 0.2
        self.alpha_multiplier = 1
        
        self.ca3_starttime = 0 #11000
        self.ca3_endtime = 612000 #21000
        self.ltp_pairings = 100
        self.ltp_freq = 100
        self.ltp_interburst_interval = 1000
        self.ltp_bursts = 1
        self.ltd_pairings = 500
        self.ltd_freq = 1
        self.stdp_pairings = 60
        self.stdp_freq = 5
        self.gaba_enabled = False
        self.plastmodifier = 1
        self.initial_weight = 1

        self.stdp_synweight = 0.5
        self.stdp_somaweight = 0.1
        self.stdp_deltat = 10

        self.epsp_measurement_enabled = True
        self.presynaptic_enabled = True
        self.full_records = True

        
    def init_object_functions(self):
        super().init_object_functions()
        h("objref vs")
        h("objref vs_ec")
        h("objref vs_gaba")
        #h("objref syntimes")
        self.default_parameters()

        
    def default_parameters(self):
        self.plasticity_model_parameters = {
            "gAMPAbar": 1.6e-4 * self.gbarmult * self.ampamult,
            "gNMDAbar": 3.0e-7 * self.gbarmult * self.nmdamult * 20.0,
            "gNR2Abar": 1 * self.nr2amult,
            "gNR2Bbar": 1 * self.nr2b_multiplier * self.nr2bmult,
            "Alpha_ampa": 1.1 * self.alpha_multiplier,
            "Alpha_nr2a": 0.5 * self.alpha_multiplier,
            "Alpha_nr2b": 0.1 * self.alpha_multiplier,
            "Beta_ampa": 0.19,
            "Beta_nr2a": 0.024,
            "Beta_nr2b": 0.0075,
            "AMPA_tau1": 0.5,
            "AMPA_tau2": 3,
            "A_ltp": 1.0e-4 * self.A_pm_mult * self.A_mult * 8.0e0,
            "A_ltd": 4.0e-3 * self.A_pm_mult * self.A_mult * 8.0e0,
            "v_thresh1": -63,
            "v_thresh2": -63,
            "v_trace1_thresh2": 0.2,
            "v_trace1_tau": 10,
            "v_trace2_tau": 10,
            "g_nmda_trace_ltp_tau": 100,
            "g_nmda_trace_ltd_tau": 3,
            "X_trace_tau": 15,
            "wmin": self.wmin,
            "wmax": self.wmax,
            "hill_midpoint_ltp": 1.7e-2,# * self.ltp_midpoint_multiplier,
            "hill_midpoint_ltp2": self.ltp_midpoint_multiplier, #9e-2,
            "hill_coef_ltp": 4,
            "hill_midpoint_ltd": 3.3e-2,
            "hill_coef_ltd": 2,
            "X_max": 1, 
            "moving_threshold_hill_ltp_multiplier": 1,
            "moving_threshold_hill_ltp_tau": 100.0,
            "moving_threshold_hill_ltp2_multiplier": 1,
            "moving_threshold_hill_ltp2_tau": 100.0,
            "moving_threshold_hill_ltd_multiplier": 1,
            "moving_threshold_hill_ltd_tau": 100.0,
        }


    def set_up_records(self, synapses, soma):
        if self.full_records:
            self.records_dendrite = { "cai": h.Vector(), }
            self.records_dendrite['cai'].record(self.monitor_dendrite._ref_cai)
            
            self.set_up_synapse_records_full(synapses, self.records_synapse)
        else:
            self.set_up_synapse_records(synapses, self.records_synapse)

        # record soma:
        self.records_soma = {"v": h.Vector()}
        self.records_soma['v'].record(soma._ref_v)
        
        # record time:
        self.records_t = h.Vector()
        self.records_t.record(h._ref_t)


    def set_up_synapse_records(self, synapses, records_synapse):
        for synapse in synapses:
            records_synapse.append(
                {
                    "i": h.Vector(),
                    "v": h.Vector(),
                    "g_nmda": h.Vector(),
                    "weight": h.Vector(),
                })
            records_synapse[-1]['v'].record(synapse.get_segment()._ref_v)
            records_synapse[-1]['i'].record(synapse._ref_i)
            records_synapse[-1]['g_nmda'].record(synapse._ref_g_nmda)
            records_synapse[-1]['weight'].record(synapse._ref_weight)


    def set_up_synapse_records_full(self, synapses, records_synapse):
        for synapse in synapses:
            records_synapse.append(
                {
                    "i": h.Vector(),
                    "v": h.Vector(),
                    "g_nmda": h.Vector(),
                    "g_nr2a": h.Vector(),
                    "g_nr2b": h.Vector(),
                    "i_nmda": h.Vector(),
                    "i_ampa": h.Vector(),
                    "v_threshed1": h.Vector(),
                    "v_threshed2": h.Vector(),
                    "g_nmda_trace_ltp": h.Vector(),
                    "g_nmda_trace_ltd": h.Vector(),
                    "w_ltp": h.Vector(),
                    "w_ltd": h.Vector(),
                    "weight": h.Vector(),
                })
            records_synapse[-1]['v'].record(synapse.get_segment()._ref_v)
            records_synapse[-1]['i'].record(synapse._ref_i)
            records_synapse[-1]['g_nmda'].record(synapse._ref_g_nmda)
            records_synapse[-1]['g_nr2a'].record(synapse._ref_g_nr2a)
            records_synapse[-1]['g_nr2b'].record(synapse._ref_g_nr2b)
            records_synapse[-1]['i_nmda'].record(synapse._ref_i_nmda)
            records_synapse[-1]['i_ampa'].record(synapse._ref_i_ampa)
            records_synapse[-1]['v_threshed1'].record(synapse._ref_v_threshed1)
            records_synapse[-1]['v_threshed2'].record(synapse._ref_v_threshed2)
            records_synapse[-1]['g_nmda_trace_ltp'].record(synapse._ref_g_nmda_trace_ltp)
            records_synapse[-1]['g_nmda_trace_ltd'].record(synapse._ref_g_nmda_trace_ltd)
            records_synapse[-1]['w_ltp'].record(synapse._ref_w_ltp)
            records_synapse[-1]['w_ltd'].record(synapse._ref_w_ltd)
            records_synapse[-1]['weight'].record(synapse._ref_weight)


    def add_plasticity_synapse(self, dendrite):        
        if os.environ["MILEDIDEBUG"] == '1':
            print("adding nmda plasticity syn")
            print(dendrite)
        synapse = self.add_synapse(dendrite, h.AMPANMDASynPlast)

        for kwarg in self.plasticity_model_parameters:
            setattr(synapse, kwarg, self.plasticity_model_parameters[kwarg])

        return synapse



    # ========================================================
    # ===================== =Analysis= =======================
    # ========================================================
    
    def get_weights(self):
        weights = []
        for item in self.records_synapse:
            weights.append(list(item["weight"])[-1])
        return weights

    
    def save_weights(self):
        weights = []
        for item in self.records_synapse:
            weights.append(list(item["weight"])[-1])
        Path(f"weights").mkdir(parents=True, exist_ok=True)

        with open(f"weights/weights.npy", 'wb') as f:
            np.save(f, weights)

            
    def mult_weights(self, multi=100):
        for i, synapse in enumerate(self.plastic_synapses):
            synapse.w = max(-0.5, synapse.w*multi)
        
        
    def load_weights(self):
        weights = np.load(f"weights/weights.npy", allow_pickle=True).tolist()
        for i, synapse in enumerate(self.plastic_synapses):
            if os.environ["MILEDIDEBUG"] == '1':
                pprint(synapse)
            synapse.w = weights[i]

