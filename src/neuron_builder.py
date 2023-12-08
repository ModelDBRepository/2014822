from .helper import neuron_node as ch
from neuron import h
from subprocess import DEVNULL, STDOUT, run, PIPE, Popen, check_call
import os
import shutil
import glob
import time
from pathlib import Path
import numpy as np
from pprint import pprint

class NeuronBuilder():

    mod_loaded = False
    root_path = os.getcwd()
    hoc_dir = os.path.join(root_path, 'nrnhoc')

    
    def __init__(self, recompile=True, silent=False, hocfile='cell_seed4_0_control_cell.hoc'):
        self.hocfile = hocfile
        self.default_parameters()
        self.silent = silent
        self.init_once(recompile)


    def default_parameters(self):
        self.neuron_parameters = {
            "syn_count": 50,
            "syn_density": 0.8,
            "rand_seed": 88,
            "somadist_L": 140,
            "trunkdist_L": None,
            "syn_diam": None,
            "gbar_kca_alz_mult": 1.0, #4.5,
            "gcalbar_cal_alz_mult": 1.0, #1.3,
            "A_mult": 1,
            "aicd_channels": 0,
            "aicd_nr2b": 0,
            "ampamult": 1,
            "nmdamult": 1,
            "nr2bmult": 1,
            "nr2bmult_aicd": 1,
            "plastmodifier": 1,
            "gaba_enabled": 0,
            "ltp_midpoint_multiplier": 1,
            "alpha_multiplier": 1,
            "syn_density_multiplier": 1,
        }
        if self.hocfile == "cell_seed4_0_AICD_cell.hoc":
            self.neuron_variable_name = "AICD_mouse_HM" 
        else:
            self.neuron_variable_name = "control_mouse_HM" 
        self.tstop = 0
        self.v_init = -70
        self.Vrest = -70
        self.celsius = 34
        self.recall_phase = True
        self.ca3_dendrites = []
        self.olm_dendrites = []
        self.basket_dendrites = []
        self.bistratified_dendrites = []
        self.ec_dendrites = []
        self.ca3_plastic_dendrites = []
        self.ec_plastic_dendrites = []
        
    def init_once(self, recompile=True):
        if NeuronBuilder.mod_loaded:
            if not self.silent:
                print("Mod files already loaded. Skipping...")
        else:
            self.compile_load_mod_files(recompile)

        hoc_path = os.path.join(self.hoc_dir, self.hocfile)
        h.load_file(str(hoc_path).replace('\\', '/'))

        h('objref cell')

        self.set_cvode(1)
        h(f'cell = new {self.neuron_variable_name}("morphology")')
        h(f"v_init={self.v_init}")
        h(f"Vrest={self.Vrest}")
        h(f"celsius={self.celsius}")

        if not self.silent:
            h('proc advance() { nrnpython("NeuronBuilder.myadvance()") }')
        h('proc init() { nrnpython("NeuronBuilder.myinit()") }')
        h('proc run() { nrnpython("NeuronBuilder.myrun()") }')
        self.set_network()


    def compile_mod_files(root_path):
        if os.environ["MILEDIDEBUG"] == '1':
            print("Recompiling neuron (running nrnivmodl)...")
        if os.name != 'nt':
            if os.path.isdir(root_path + '/mechanisms/x86_64'):
                shutil.rmtree(root_path + '/mechanisms/x86_64')
        else:
            for filePath in glob.glob(root_path + "/mechanisms/*.dll", recursive=False):
                try:
                    os.remove(filePath)
                except:
                    print("Error while deleting file : ", filePath)
        os.chdir(root_path + '/mechanisms')
        compile_program = "nrnivmodl"
        if os.name == 'nt':
            compile_program = "nrnivmodl.bat"
        p = check_call(compile_program)
        os.chdir(root_path)
        if os.environ["MILEDIDEBUG"] == '1':
            print("Compiled neuron. Loading dll...")


    def load_mod_files(root_path):
        libpath = "/mechanisms/x86_64/.libs/libnrnmech.so"
        if os.name == 'nt':
            libpath = "/mechanisms/nrnmech.dll"

        h.nrn_load_dll(str(root_path + libpath))

        if os.environ["MILEDIDEBUG"] == '1':
            print("Loaded dll.")


    def compile_load_mod_files(self, recompile=True):
        if (recompile):
            NeuronBuilder.compile_mod_files(self.root_path)
        else:
            if os.environ["MILEDIDEBUG"] == '1':
                print("Skipping compiling neuron. Loading dll...")
            
        NeuronBuilder.load_mod_files(self.root_path)
        NeuronBuilder.mod_loaded = True


    def set_alzheimers(self, aicd_channels=0, aicd_nr2b=0):
        self.neuron_parameters["aicd_channels"] = aicd_channels
        self.neuron_parameters["nr2bmult_aicd"] = aicd_nr2b

        
    def set_cvode(self, val):
        h(f"_ = cvode_active({val})")
        h.using_cvode_ = 1 if val else 0
        if val:
            h.CVode()
            h.CVode().active(True)

    def set_parameters(self, kwargs):
        if os.environ["MILEDIDEBUG"] == '1':
            print("Setting neuron parameters")
        self.default_parameters()
        for kwarg in kwargs:

            if os.environ["MILEDIDEBUG"] == '1':
                print(f"- [{kwarg}]: {kwargs[kwarg]}")
            self.neuron_parameters[kwarg] = kwargs[kwarg]


    def set_tstop(self, tstop):
        self.tstop = tstop
        h.tstop = tstop

        
    def run(self):
        h.tstop = self.tstop
        if os.environ["MILEDIDEBUG"] == '1':
            print(f"Running neuron simulation till {h.tstop}")

        if os.environ["MILEDIDEBUG"] == '1':
            self.spikebuilder.pprintinfo()
        h.init()
        h.run()
        
        if not self.silent:
            NeuronBuilder.printtime()


    @staticmethod
    def myrun():
        h.running_ = 1
        h("debug_counter = 0")
        h(f"start_time = {time.time()}")

        h.stdinit()
        h.cvode_simgraph()
        h.realtime = 0
        _ = h.setdt()
        h.init()

        h.continuerun(h.tstop)

        
    @staticmethod
    def myadvance():
        h.fadvance()
        h("debug_counter += 1")
        if h.debug_counter % 1000 == 0:
            NeuronBuilder.printtime()
        
    @staticmethod  
    def printtime():
        ratio_done = h.t/h.tstop
        ratio_remaining = 1 - ratio_done
        time_elapsed = time.time() - h.start_time
        time_for_1_percent = time_elapsed/ratio_done
        time_left = int(time_for_1_percent * ratio_remaining)
        hours = int(time_left / 3600)
        minutes = int((time_left-(hours*3600)) / 60)
        seconds = int(time_left % 60)
        print(f"{h.debug_counter}: {h.t:.2f}/{h.tstop:.2f} ({(100*ratio_done):.2f}%). Time remaining: {hours:02d}:{minutes:02d}:{seconds:02d}")

        
    @staticmethod
    def myinit():
        h.t = 0
        h.finitialize(h.Vrest)
        if h.cvode.active():
            h.cvode.re_init()
        else:
            h.fcurrent()

        # TODO
        # reset cai and other variables here to init
        h.frecord_init()


    # ===================================================================
    # ===================================================================
    # ===================================================================
    
    def set_network(self):

        if os.environ["MILEDIDEBUG"] == '1':
            print("Building morphology network...")
        self.soma = h.cell.soma[0]
        self.all_nodes = ch.buildNodes(self.soma)
        self.root_node_object = self.all_nodes[str(self.soma(0.5))]

        if os.environ["MILEDIDEBUG"] == '1':
            print("Building morphology network done.")

        
    def load_session(self, sesfile="main.ses"):
        h.load_file(f"{sesfile}")

        
    def aicd_mult(self, multi):
        return 1 + (multi - 1) * self.neuron_parameters["aicd_channels"]

    
    def add_coupling(self, gbar_kca):
        coupling_multiplier = 1

        for dend in self.dendrite_cluster:
            dend.item.gbar_kca = gbar_kca * coupling_multiplier

        for dend in self.ca3_dendrites:
            dend[0].gbar_kca = gbar_kca * coupling_multiplier

        for dend in self.ec_dendrites:
            dend[0].gbar_kca = gbar_kca * coupling_multiplier

        for dend in self.ca3_plastic_dendrites:
            dend[0].item.gbar_kca = gbar_kca * coupling_multiplier

        for dend in self.ec_plastic_dendrites:
            dend[0].item.gbar_kca = gbar_kca * coupling_multiplier
            
    def adjust_dendrite_mechanisms(self):
        gcalbar_cal_alz_multiplier = self.aicd_mult(self.neuron_parameters["gcalbar_cal_alz_mult"])
        gbar_kca_multiplier = self.aicd_mult(self.neuron_parameters["gbar_kca_alz_mult"])

        if os.environ["MILEDIDEBUG"] == '1':
            print(f"gcalbar_cal_alz_mult: {gcalbar_cal_alz_multiplier}")
            print(f"gbar_kca_multiplier: {gbar_kca_multiplier}")

        if self.hocfile == "cell_seed4_0_AICD_cell.hoc":
            gcalbar_cal = 7.7381721218383238e-07 * gcalbar_cal_alz_multiplier
            gbar_kca = 8.1047396376148008e-05 * gbar_kca_multiplier
        else:
            gcalbar_cal = 4.1962096014411989e-07 * gcalbar_cal_alz_multiplier
            gbar_kca = 8.1047396376148008e-05 * gbar_kca_multiplier


        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'gcalbar_cal'):
                nobj.gcalbar_cal = gcalbar_cal
            if hasattr(nobj, 'gbar_kca'):
                nobj.gbar_kca = gbar_kca

        self.add_coupling(gbar_kca)
        
    def zero_ih(self):
        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'ghdbar_hd'):
                nobj.ghdbar_hd = 0
            if hasattr(nobj, 'e_pas'):
                nobj.e_pas = -70

                
    def zero_kad(self):
        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'gkabar_kad'):
                nobj.gkabar_kad = nobj.gkabar_kad*0.8


    def set_beta_amyloid(self, beta, beta_midpoint):
        ampa = 1.0
        # alpha = 1 + (1.2 - 1) * beta
        syn_density = 1 + (0.8 - 1) * beta

        #ampa = 1 + (0.75 - 1) * beta # 0.75
        #ampa = 1 + (0.6 - 1) * beta # 0.6
        #alpha = 1 + (4.0 - 1) * beta # 2.0
        #syn_density = 1 + (0.75 - 1) * beta # 0.75

        if beta > 0:
            alpha = 1.2
        else:
            alpha = 1.0

        self.neuron_parameters["ampamult"] = ampa
        self.neuron_parameters["ltp_midpoint_multiplier"] = beta_midpoint # 1.75
        self.neuron_parameters["alpha_multiplier"] = alpha
        self.neuron_parameters["syn_density_multiplier"] = syn_density

    def set_ih(self, val):
        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'ghdbar_hd'):
                nobj.ghdbar_hd = val *  nobj.ghdbar_hd

                
    def set_na(self, val):
        for node in self.all_nodes:
            if "axon" in node:
                continue
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'gbar_nax'):
                nobj.gbar_nax = val *  nobj.gbar_nax

                
    def set_kdrbar(self, val):
        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'gkdrbar_kdr'):
                nobj.gkdrbar_kdr = val *  nobj.gkdrbar_kdr

                
    def set_kadbar(self, val):
        for node in self.all_nodes:
            nobj = self.all_nodes[node].item
            if hasattr(nobj, 'gkabar_kad'):
                nobj.gkabar_kad = val * nobj.gkabar_kad
                
                
    def printe(self):
        h('forall print secname(), e_pas')
    
    # ========================================================
    # ==================== =Protocols= =======================
    # ========================================================

    
    def build_protocol(self, spikebuilder, protocol_function, **kwargs):
        self.spikebuilder = spikebuilder
        spikebuilderparameters = {
                "ampamult": self.neuron_parameters["ampamult"],
                "nr2bmult": self.neuron_parameters["nr2bmult"] * self.neuron_parameters["nr2bmult_aicd"],
                "plastmodifier": self.neuron_parameters["plastmodifier"],
                "gaba_enabled": self.neuron_parameters["gaba_enabled"],
                "nmdamult": self.neuron_parameters["nmdamult"],
                "A_mult": self.neuron_parameters["A_mult"],
                "ltp_midpoint_multiplier": self.neuron_parameters["ltp_midpoint_multiplier"],
                "alpha_multiplier": self.neuron_parameters["alpha_multiplier"],
                }
        if os.environ["MILEDIDEBUG"] == '1':
            print("Building spiking protocol. Spikebuilder parameters:")
            print(f'nr2bmult: {self.neuron_parameters["nr2bmult"]}, nr2bmult_aicd: {self.neuron_parameters["nr2bmult_aicd"]}')
            print(spikebuilderparameters)
            print(self.neuron_variable_name)

        spikebuilder.set_parameters(spikebuilderparameters)
        spikebuilder.set_alzheimers(self.neuron_parameters["aicd_nr2b"])

        syn_count = self.neuron_parameters["syn_count"]
        syn_density = self.neuron_parameters["syn_density"] * self.neuron_parameters["syn_density_multiplier"]

        if os.environ["MILEDIDEBUG"] == '1':
            print(f"Synapse count: {syn_count}")
        
        if self.neuron_parameters["syn_count"] > 0:
            (dendrite_cluster, centernode) = ch.create_random_cluster(self.all_nodes, syn_density, syn_count, self.neuron_parameters["somadist_L"], self.neuron_parameters["trunkdist_L"], self.neuron_parameters["syn_diam"], self.neuron_parameters["rand_seed"])
            monitor_dend = centernode.item
            self.dendrite_cluster = dendrite_cluster
            #pprint(dendrite_cluster)
        else:
            dendrite_cluster = []
            self.dendrite_cluster = []
            monitor_dend = None
        
        h('objref monitor_synapse')
        h('objref monitor_dendrite')
        h('objref somagaba')
        if monitor_dend:
            spikebuilder.set_monitor_dendrite(monitor_dend)
            h.monitor_dendrite = monitor_dend

        totalruntime = protocol_function(dendrite_cluster, self.soma(0.5), **kwargs)
        
        self.tstop = totalruntime
        h.tstop = self.tstop

        # if "AICD" in self.neuron_variable_name:
        #     self.adjust_dendrite_mechanisms()


        if len(spikebuilder.plastic_synapses) > 0:
            h.monitor_synapse = spikebuilder.plastic_synapses[0]

        if len(spikebuilder.gaba_synapses) > 0:
            h.somagaba = spikebuilder.gaba_synapses[0]

        if os.environ["MILEDIDEBUG"] == '1':
            print(f"Building spiking protocol done. tstop: {self.tstop}.")

        
    # ========================================================
    # ======================= =Misc= =========================
    # ========================================================

    
    def add_iclamp(self, spikebuilder, current, nodes=[], clampall=False):
        self.set_cvode(0)
        
        chi = self.root_node_object.children[-1]
        for i in range(20):
            chi = chi.children[-1]
            print(chi)



        spikebuilder.set_alzheimers(self.neuron_parameters["aicd_nr2b"])
        if clampall:
            spikebuilder.add_iclamp(self.all_nodes.values(), self.soma(0.5), current)
        else:
            spikebuilder.add_iclamp(nodes, chi.item, current)
            #spikebuilder.add_iclamp(nodes, self.soma(0.5), current)

            
    def add_vclamp(self, spikebuilder, voltage, nodes=[], clampall=False):
        self.set_cvode(0)
        spikebuilder.set_alzheimers(self.neuron_parameters["aicd_nr2b"])
        if clampall:
            spikebuilder.add_vclamp(self.all_nodes.values(), self.soma(0.5), voltage)
        else:
            spikebuilder.add_vclamp(nodes, self.soma(0.5), voltage)
        h('objref soma_vclamp')
        h.soma_vclamp = spikebuilder.clamp_connections[-1]


    def set_ec(self, node_names):
        for (node_name, w) in node_names:
            self.ec_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_ec_dend')
            h.monitor_ec_dend = self.ec_dendrites[0][0].item

    def set_olm(self, node_names):
        for (node_name, w) in node_names:
            self.olm_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_olm_dend')
            h.monitor_olm_dend = self.olm_dendrites[0][0].item

    def set_bistratified(self, node_names):
        for (node_name, w) in node_names:
            self.bistratified_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_bistratified_dend')
            h.monitor_bistratified_dend = self.bistratified_dendrites[0][0].item
            
    def set_basket(self, node_names):
        for (node_name, w) in node_names:
            self.basket_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_basket_dend')
            h.monitor_basket_dend = self.basket_dendrites[0][0].item

    def set_ca3(self, node_names):
        for (node_name, w) in node_names:
            self.ca3_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_ca3_dend')
            h.monitor_ca3_dend = self.ca3_dendrites[0][0].item
    
    def set_ca3_plastic(self, node_names):
        for (node_name, w) in node_names:
            self.ca3_plastic_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_ca3_plastic_dend')
            h.monitor_ca3_plastic_dend = self.ca3_plastic_dendrites[0][0].item

    def set_ec_plastic(self, node_names):
        for (node_name, w) in node_names:
            self.ec_plastic_dendrites.append((self.all_nodes[node_name], w))
        
        if node_names:
            h('objref monitor_ec_plastic_dend')
            h.monitor_ec_plastic_dend = self.ec_plastic_dendrites[0][0].item

    # ========================================================
    # ===================== =External= =======================
    # ========================================================

    
    def build_protocol_ltp(self, spikebuilder):
        self.build_protocol(spikebuilder, spikebuilder.build_protocol_ltp)

        
    def build_protocol_ltd(self, spikebuilder):
        self.build_protocol(spikebuilder, spikebuilder.build_protocol_ltd)

        
    def build_protocol_ca3(self, spikebuilder):
        self.build_protocol(spikebuilder, spikebuilder.build_protocol_ca3, ca3id=1, ec_dendrite=self.ec_dend)


    def build_protocol_theta_gamma(self, spikebuilder):
        self.neuron_parameters["syn_count"] = 0


        # CA3 plastic
        self.set_ca3_plastic(
            [
                (f"{self.neuron_variable_name}[0].apic[21](0.5)", 1.0),
            ]
        )

        # CA3 non plastic 
        self.set_ca3(
            [
                (f'{self.neuron_variable_name}[0].apic[21](0.5)', 1.0),
            ]
        )

        # EC plastic
        self.set_ec_plastic(
            [
                #f"{self.neuron_variable_name}[0].apic[39](0.5)",
            ]
        )

        # EC non plastic
        self.set_ec(
            [
                (f'{self.neuron_variable_name}[0].apic[41](0.5)', 1.0),
            ]
        )
       
        # Basket
        self.set_basket(
            [
                (f'{self.neuron_variable_name}[0].soma[0](0.5)', 1.0),
            ]
        )
        
        # Bistratified
        self.set_bistratified(
            [
                (f'{self.neuron_variable_name}[0].apic[7](0.5)', 1.0),
            ]
        )
        
        # OLM
        self.set_olm(
            [
                (f'{self.neuron_variable_name}[0].apic[38](0.5)', 1.0),
            ]
        )


        self.build_protocol(spikebuilder, spikebuilder.build_protocol_theta_gamma, 
            recall_phase=self.recall_phase,
            ec_dendrites=self.ec_dendrites,
            olm_dendrites=self.olm_dendrites,
            bistratified_dendrites=self.bistratified_dendrites,
            basket_dendrites=self.basket_dendrites,
            ca3_dendrites=self.ca3_dendrites,
            ca3_plastic_dendrites=self.ca3_plastic_dendrites,
            ec_plastic_dendrites=self.ec_plastic_dendrites,
        )

    def build_protocol_stdp(self, spikebuilder):
        self.build_protocol(spikebuilder, spikebuilder.build_protocol_stdp)


    def build_protocol_ratio(self, spikebuilder):
        self.build_protocol(spikebuilder, spikebuilder.build_protocol_ratio)

    # ========================================================
    # ====================== =Debug= =========================
    # ========================================================

    def write_debug_files(self, folder):
        Path(f"debug/{folder}").mkdir(parents=True, exist_ok=True)

        with open(f"debug/{folder}/neuron_parameters.npy", 'wb') as f:
            np.save(f, self.neuron_parameters)
