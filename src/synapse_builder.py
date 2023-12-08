from neuron import h
from pprint import pprint
from pathlib import Path
import numpy as np
import pickle
import time
import pathlib
import os

class SynapseBuilder():

    epsp_start = 200
    epsp_end = 2000
    after_epsp_wait = 200
    stimulation_start = 500
    
    root_path = os.getcwd()
    
    def __init__(self):
        self.init_object_vars()
        self.init_object_functions()
        
    def init_object_vars(self):
        self.democracy = False
        self.records_synapse = []
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
        self.gamma_pairings = 5
        self.theta_pairings = 20
        self.gamma_freq = 80
        self.theta_freq = 4
        self.gbarmult = 1
        self.nr2b_cure = 1
        self.A_pm_mult = 10
        self.A_mult = 1
        self.gaba_weight = 1.0
        self.ampamult = 1
        self.nmdamult = 1
        self.nr2amult = 1
        self.nr2bmult = 1
        self.wmax = 0.5
        self.wmin = -0.5
        
        self.ca3_starttime = 0 #11000
        self.ca3_endtime = 612000 #21000
        self.ltp_pairings = 100
        self.ltp_freq = 100
        self.ltp_interburst_interval = 1000
        self.ltp_bursts = 1
        self.ltd_pairings = 100
        self.ltd_freq = 1
        self.stdp_pairings = 30
        self.stdp_freq = 5
        self.gaba_enabled = False
        self.plastmodifier = 1
        self.initial_weight = 1
        self.ltp_midpoint_multiplier = 1
        
        self.stdp_synweight = 0.5
        self.stdp_somaweight = 0.1
        self.stdp_deltat = 10

        self.epsp_measurement_enabled = False
        self.presynaptic_enabled = True
        
    def __del__(self): 
        self.reset_object_vars()
        
    def reset_object_vars(self):
        for item in self.stimulations:
            item = None
        for item in self.clamp_connections:
            item = None
        for item in self.network_connections:
            item = None
        for item in self.plastic_synapses:
            item = None
        for item in self.nonplastic_synapses:
            item = None
        for item in self.gaba_synapses:
            item = None
        for item in self.artificial_cells:
            item = None
            
        self.clamp_connections = []
        self.stimulations = []
        self.network_connections = []
        self.plastic_synapses = []
        self.nonplastic_synapses = []
        self.artificial_cells = []
        self.gaba_synapses = []
        
    def init_object_functions(self):
        h("objref vs")
        self.default_parameters()

        
    def set_parameters(self, kwargs):
        if os.environ["MILEDIDEBUG"] == '1':
            print("Setting synapse parameters")
        for kwarg in kwargs:
            if os.environ["MILEDIDEBUG"] == '1':
                print(f"- [{kwarg}]: {kwargs[kwarg]}")
            setattr(self, kwarg, kwargs[kwarg])
        self.default_parameters()

            
    def set_plasticity_model_parameters(self, kwargs):
        self.default_parameters()
        for kwarg in [x for x in kwargs]:
            self.plasticity_model_parameters[kwarg] = kwargs[kwarg]


    def default_parameters(self):
        self.plasticity_model_parameters = {
            "gAMPAbar": 1.27e-5 * self.gbarmult * self.ampamult,
            "gNMDAbar": 1.0e-5 * self.gbarmult * self.nmdamult,
            "gNR2Abar": 1 * self.nr2amult,
            "gNR2Bbar": 1 * self.nr2b_multiplier * self.nr2bmult,
            "Alpha_nr2a": 0.5,
            "Alpha_nr2b": 0.1,
            "Beta_nr2a": 0.024,
            "Beta_nr2b": 0.0075,
            "tau1": 0.1,
            "tau2": 3,
            "A_p": 4e-5 * self.A_pm_mult * self.A_mult,
            "A_m": 5e-5 * self.A_pm_mult * self.A_mult,
            "tetam": -69,
            "tetap": -64,
            "tau_0": 10,
            "tau_y": 110,
            "nmt_tau": 10,
            "X_a": -5,
            "X_b": 5000,
            "nmda_thresh_p": 1e-7,
            "nmda_thresh_d": 1e-7,
            "gnmda_tanh_b": 1e4,
            "wmin": self.wmin,
            "wmax": self.wmax,
            "gbar": 0,
            "gbarmid": (self.initial_weight - 1)/18,
            "plastmodifier": self.plastmodifier,
        }
        self.plasticity_model_parameters = {}

        
    def aicd_mult(multi, aicd):
        return 1 + (multi-1) * aicd

    
    def set_alzheimers(self, aicd = 0):
        self.nr2b_multiplier = SynapseBuilder.aicd_mult(4, aicd)

        if os.environ["MILEDIDEBUG"] == '1':
            print(f"NR2B Multiplier: {self.nr2b_multiplier}")
        self.default_parameters()
            
            
    # ========================================================
    # ===================== =Builders= =======================
    # ========================================================
    
        
    def add_network_connection(self, stimulation, synapse, w):
        if (os.environ.get("NRN_DEBUG") == "1"):
            print(f"Adding connection: {synapse}")
        netcon = h.NetCon(stimulation, synapse, 0, 0, w)
        self.network_connections.append(netcon)

                
    def add_stimulation(self, pairings, freq, start, noise=0):
        stimulation = h.NetStim()
        stimulation.number = pairings
        stimulation.interval = 1000.0 / freq
        stimulation.start = start
        stimulation.noise = noise

        if (os.environ.get("NRN_DEBUG") == "1"):
            print(f"Adding stimulation: {pairings} pair * {stimulation.interval} interval ({freq} freq). Start at: {start}")

        self.stimulations.append(stimulation)
        return stimulation

    
    # ========================================================
    # ===================== =Synapses= =======================
    # ========================================================

    
    def add_nonplastic_synapse(self, dendrite, plasticity_model):
        synapse = plasticity_model(dendrite)
        self.nonplastic_synapses.append(synapse)
        return synapse

    
    def add_synapse(self, dendrite, plasticity_model):
        synapse = plasticity_model(dendrite)
        self.plastic_synapses.append(synapse)
        return synapse

    
    def add_gaba_synapse(self, dendrite):
        synapse = self.add_nonplastic_synapse(dendrite, lambda x: h.Exp2Syn(0.5, x))
        self.gaba_synapses.append(synapse)
        synapse.tau1 = 35
        synapse.tau2 = 100
        synapse.e = -75

        return synapse


    def add_ampa_exp2syn_synapse(self, dendrite):
        synapse = self.add_nonplastic_synapse(dendrite, lambda x: h.Exp2Syn(0.5, x))
        synapse.tau1 = 0.1
        synapse.tau2 = 3
        synapse.e = 0

        return synapse

    
    def add_artificial_gaba(self, stim_in, somagaba, w):
        artificial_cell = h.IntFire4()
        self.artificial_cells.append(artificial_cell)
        self.add_network_connection(stim_in, artificial_cell, 0.4)
        for somagaba_syn in somagaba:
            self.add_network_connection(artificial_cell, somagaba_syn, w)
        

    def add_plasticity_synapse(self, dendrite):        
        if os.environ["MILEDIDEBUG"] == '1':
            print("adding nmda syn")
        synapse = self.add_synapse(dendrite, h.NMDASyn)

        for kwarg in self.plasticity_model_parameters:
            setattr(synapse, kwarg, self.plasticity_model_parameters[kwarg])

        return synapse

    
    def add_AMPA_democracy(self, synapses, weights):
        if not self.democracy:
            return
        for (i, synapse) in enumerate(synapses):
            synapse.gAMPAbar = synapse.gAMPAbar * weights[i]
    
    
    def build_synapses(self, dendrites, synapse_function):
        synapses = []
        for dendrite in dendrites:
            synapse = synapse_function(dendrite)
            synapses.append(synapse)
        return synapses

    
    def connect_iclamp(self, dendrites, current):
        for dendrite in dendrites:
            self.add_iclamp_connection(dendrite, current)

            
    def connect_vclamp(self, dendrites, voltage):
        for dendrite in dendrites:
            self.add_vclamp_connection(dendrite, voltage)

            
    def connect_synapses_stimulation(self, synapses, stimulation, w=1):
        if type(w) is list:
            for i, synapse in enumerate(synapses):
                self.add_network_connection(stimulation, synapse, w[i] * self.weight_factor)
        else:
            for synapse in synapses:
                self.add_network_connection(stimulation, synapse, w * self.weight_factor)
                
    def add_gaba_synapses(self, dendrite, count):
        synapses = []
        for i in range(count):
            synapse = self.add_gaba_synapse(dendrite)
            synapses.append(synapse)
        return synapses

            
    def build_plasticity_synapses(self, dendrites):
        return self.build_synapses(dendrites, self.add_plasticity_synapse)
    
    
    # ========================================================
    # =================== =Stimulations= =====================
    # ========================================================

    
    def build_stimulations_ltp_gaba(self):
        stimulations = []
        freq = 100
        pairings = 100
        start = self.stimulation_start
        stimulation = self.add_stimulation(pairings, freq, start-1)
        stimulations.append(stimulation)
        stimulation = self.add_stimulation(pairings, freq, start+1)
        stimulations.append(stimulation)
        stimulation_time = (1000/freq) * pairings
        interburst_time = 1000
        
        start = start + stimulation_time + interburst_time
        
        stimulation = self.add_stimulation(pairings, freq, start-1)
        stimulations.append(stimulation)
        stimulation = self.add_stimulation(pairings, freq, start+1)
        stimulations.append(stimulation)

        return stimulations


    def build_stimulations_ratio(self):
        stimulations = []
        freq = 100
        pairings = 1
        start = self.stimulation_start
        stimulation = self.add_stimulation(pairings, freq, start)
        stimulations.append(stimulation)
        stimulation_time = (1000/freq) * pairings
            
        return (stimulations, stimulation_time)


    def build_stimulations_ltp(self):
        stimulations = []
        freq = self.ltp_freq
        pairings = self.ltp_pairings
        start = self.stimulation_start
        stimulation = self.add_stimulation(pairings, freq, start)
        stimulations.append(stimulation)
        stimulation_time = (1000/freq) * pairings
    
        for i in range(self.ltp_bursts-1):

            if os.environ["MILEDIDEBUG"] == '1':
                print(f"building burst: {i}")
            start = start + stimulation_time + self.ltp_interburst_interval
            stimulation = self.add_stimulation(pairings, freq, start)
            stimulations.append(stimulation)

        stimulation_time = stimulation_time * self.ltp_bursts + self.ltp_interburst_interval * (self.ltp_bursts - 1)
        return (stimulations, stimulation_time)


    def build_stimulations_ltd(self):
        stimulations = []
        freq = self.ltd_freq
        pairings = self.ltd_pairings
        start = self.stimulation_start
        stimulation = self.add_stimulation(pairings, freq, start)
        stimulations.append(stimulation)
        stimulation_time = (1000/freq) * pairings

        return (stimulations, stimulation_time)

    
    def build_epsp_measurement(self, synapses, totalruntime):
        stimulations = self.build_stimulations_epsp(totalruntime)

        for stimulation in stimulations:
            self.connect_synapses_stimulation(synapses, stimulation)

        return stimulations
            
    def build_stimulations_epsp(self, totalruntime):
        stimulations = []
        freq = 1
        pairings = 1
        start = self.epsp_start
        stimulation = self.add_stimulation(pairings, freq, start)
        stimulations.append(stimulation)

        start = totalruntime - self.after_epsp_wait
        stimulation = self.add_stimulation(pairings, freq, start)
        stimulations.append(stimulation)

        return stimulations

            
    # ========================================================
    # ======================= =Misc= =========================
    # ========================================================

    def set_up_records(self, synapses, soma):
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
                })
            records_synapse[-1]['v'].record(synapse.get_segment()._ref_v)
            records_synapse[-1]['i'].record(synapse._ref_i)


    def get_parameter(self, synapseid, parameter, t):
        tid = next(x[0] for x in enumerate(self.records_t) if x[1] > t)

        rec = self.records_synapse[synapseid][parameter][tid]
        return (self.records_t[tid], rec)

    
    def get_parameter_range(self, synapseid, parameter, tfrom, tto):
        tid_from = next(x[0] for x in enumerate(self.records_t) if x[1] > tfrom)
        tid_to = next(x[0] for x in enumerate(self.records_t) if x[1] > tto)

        data = list(self.records_synapse[synapseid][parameter])
        rec = data[tid_from:tid_to]
        return (self.records_t[tid_from], self.records_t[tid_to], rec)

        
    def add_vclamp(self, all_nodes, soma, voltage):
        dendrites = [x.item for x in all_nodes]
        
        self.connect_vclamp(dendrites, voltage)
        self.connect_vclamp([soma], voltage)


    def add_iclamp(self, all_nodes, soma, current):
        dendrites = [x.item for x in all_nodes]
        
        self.connect_iclamp(dendrites, current)
        self.connect_iclamp([soma], current)

        
    def set_monitor_dendrite(self, monitor_dendrite):
        self.monitor_dendrite = monitor_dendrite

        
    # ========================================================
    # ==================== =Protocols= =======================
    # ========================================================

    def build_protocol_ratio(self, dendrite_nodes, soma):
        dendrites = [x.item for x in dendrite_nodes]
        synapses = self.build_plasticity_synapses(dendrites)

        self.set_up_records(synapses, soma)
        self.add_AMPA_democracy(synapses, [x.weight_by_distance for x in dendrite_nodes])
        
        (stimulations, stimulation_time) = self.build_stimulations_ratio()
        for stimulation in stimulations:
            self.connect_synapses_stimulation(synapses, stimulation)          

        self.totalruntime = stimulation_time + self.stimulation_start + self.after_epsp_wait

        return self.totalruntime

    
    def build_protocol_ltp(self, dendrite_nodes, soma):
        dendrites = [x.item for x in dendrite_nodes]
        synapses = self.build_plasticity_synapses(dendrites)
        gabasynapses = self.add_gaba_synapses(soma, 1)

        self.set_up_records(synapses, soma)
        self.add_AMPA_democracy(synapses, [x.weight_by_distance for x in dendrite_nodes])
        
        (stimulations, stimulation_time) = self.build_stimulations_ltp()
        for stimulation in stimulations:
            self.connect_synapses_stimulation(synapses, stimulation)
            if self.gaba_enabled:
                self.connect_synapses_stimulation(gabasynapses, stimulation, self.gaba_weight)
                #self.add_artificial_gaba(stimulation, gabasynapses, self.gaba_weight)
                


        if self.epsp_measurement_enabled:

            if os.environ["MILEDIDEBUG"] == '1':
                print(f"stimulation_time: {stimulation_time}")
            totalruntime = stimulation_time + self.stimulation_start + self.epsp_end + self.after_epsp_wait

            if os.environ["MILEDIDEBUG"] == '1':
                print(f"Measurement TotalTime: {totalruntime}")
            self.build_epsp_measurement(synapses, totalruntime)
        else:
            totalruntime = stimulation_time + self.stimulation_start + self.after_epsp_wait

        self.totalruntime = totalruntime

        return totalruntime

    
    def build_protocol_ltd(self, dendrite_nodes, soma):
        dendrites = [x.item for x in dendrite_nodes]
        synapses = self.build_plasticity_synapses(dendrites)

        self.set_up_records(synapses, soma)
        self.add_AMPA_democracy(synapses, [x.weight_by_distance for x in dendrite_nodes])
        
        (stimulations, stimulation_time) = self.build_stimulations_ltd()
        for stimulation in stimulations:
            self.connect_synapses_stimulation(synapses, stimulation)

        if self.epsp_measurement_enabled:
            if os.environ["MILEDIDEBUG"] == '1':
                print(f"stimulation_time: {stimulation_time}")
            totalruntime = stimulation_time + self.stimulation_start + self.epsp_end + self.after_epsp_wait
            if os.environ["MILEDIDEBUG"] == '1':
                print(f"Measurement TotalTime: {totalruntime}")
            self.build_epsp_measurement(synapses, totalruntime)
        else:
            totalruntime = stimulation_time + self.stimulation_start + self.after_epsp_wait

        self.totalruntime = totalruntime
        return totalruntime

    # ========================================================
    # ===================== =Analysis= =======================
    # ========================================================
    
    def get_weights(self):
        weights = []
        for item in self.records_synapse:
            weights.append(list(item["synweight"])[-1])
        return weights

    
    def save_weights(self):
        weights = []
        for item in self.records_synapse:
            weights.append(list(item["gbar"])[-1])
        Path(f"weights").mkdir(parents=True, exist_ok=True)

        with open(f"weights/weights.npy", 'wb') as f:
            np.save(f, weights)

            
    def mult_weights(self, multi=100):
        for i, synapse in enumerate(self.plastic_synapses):
            synapse.gbar = max(-0.5, synapse.gbar*multi)
        
        
    def load_weights(self):
        weights = np.load(f"weights/weights.npy", allow_pickle=True).tolist()
        for i, synapse in enumerate(self.plastic_synapses):
            pprint(synapse)
            synapse.gbar = weights[i]


    def print_epsp_change(self):
        epspchange = self.get_epsp_change()
        if epspchange is not None:
            print(f"EPSP change: {epspchange:6.2f}%")

    
    def get_epsp_change(self):
        if not self.epsp_measurement_enabled:
            print("EPSP measurement is not enabled.")
            return None

        if (h.t < int(self.totalruntime - self.after_epsp_wait + 50)):
            print("Stopped early. Can't calculate EPSP.")
            return None

        min_start1 = sum([x <= (self.epsp_start - 50) for x in list(self.records_t)])
        min_end1 = sum([x <= self.epsp_start for x in list(self.records_t)])

        max_start1 = sum([x <= self.epsp_start for x in list(self.records_t)])
        max_end1 = sum([x <= (self.epsp_start + 50) for x in list(self.records_t)])

        min_start2 = sum([x <= (self.totalruntime - self.after_epsp_wait - 50) for x in list(self.records_t)])
        min_end2 = sum([x <= (self.totalruntime - self.after_epsp_wait) for x in list(self.records_t)])

        max_start2 = sum([x <= (self.totalruntime - self.after_epsp_wait) for x in list(self.records_t)])
        max_end2 = sum([x <= (self.totalruntime - self.after_epsp_wait + 50) for x in list(self.records_t)])

        min1 = min(list(self.records_soma["v"])[min_start1:min_end1])
        max1 = max(list(self.records_soma["v"])[max_start1:max_end1])

        min2 = min(list(self.records_soma["v"])[min_start2:min_end2])
        max2 = max(list(self.records_soma["v"])[max_start2:max_end2])
        
        
        before = (max1 - min1)
        after = (max2 - min2)
        change = (after * 100) / before

            
        return change


    def write_debug_files(self, folder="0"):
        Path(f"debug/{folder}").mkdir(parents=True, exist_ok=True)

        if self.records_ec_synapse is not None:
            with open(f"debug/{folder}/ec_synapse.npy", 'wb') as f:
                np.save(f, self.records_ec_synapse)

        with open(f"debug/{folder}/synapse.npy", 'wb') as f:
            np.save(f, self.records_synapse)

        with open(f"debug/{folder}/dendrite.npy", 'wb') as f:
            np.save(f, self.records_dendrite)
            
        with open(f"debug/{folder}/soma.npy", 'wb') as f:
            np.save(f, self.records_soma)
            
        with open(f"debug/{folder}/t.npy", 'wb') as f:
            np.save(f, self.records_t)

        with open(f"debug/{folder}/epsp.npy", 'wb') as f:
            np.save(f, self.get_epsp_change())

        with open(f"debug/{folder}/parameters.npy", 'wb') as f:
            np.save(f, self.plasticity_model_parameters)

        with open(f"debug/{folder}/spikecount.npy", 'wb') as f:
            np.save(f, self.calculate_spikes())

    def calculate_spikes_soma(self):
        spikestarted = False
        spikecount = 0
        for vmem in self.records_soma['v']:
            if vmem > -30 and not spikestarted:
                spikestarted = True
                spikecount += 1
            if vmem < -30 and spikestarted:
                spikestarted = False
        return spikecount

    def calculate_spikes_synapse(self):
        spikes = []
        for syndata in self.records_synapse:
            spikestarted = False
            spikecount = 0
            for vmem in syndata['v']:
                if vmem > -30 and not spikestarted:
                    spikestarted = True
                    spikecount += 1
                if vmem < -30 and spikestarted:
                    spikestarted = False
            spikes.append(spikecount)
        return spikes

    def calculate_spikes(self):
        synapse_spikes = self.calculate_spikes_synapse()
        soma_spikes = self.calculate_spikes_soma()

        return {"synapse_spikes": synapse_spikes, "soma_spikes": soma_spikes }
    

    def get_parameters(self):
        return self.plasticity_model_parameters
        

    def pprintinfo(self):
        pprint(self.get_parameters())
