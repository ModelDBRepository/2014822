#!/usr/bin/python
from src.synapse_builder import SynapseBuilder
from src.neuron_builder import NeuronBuilder
from src.synapse_builder_ampanmdasynplast import AMPANMDASynPlastBuilder
from src.helper.config import load_config, save_config
from multiprocessing import Process
from datetime import datetime
import argparse
import os
import numpy as np
from glob import glob

SynapseType = AMPANMDASynPlastBuilder

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)

def parse_args():
    parser = argparse.ArgumentParser(description='Build and run neuron simulation',
                                    formatter_class=SmartFormatter)

    parser.add_argument("-p", "--protocol",
                            choices=[
                                'all',
                                ],
                            default="all",
                            help="Protocol to run"
                            )
    
    parser.add_argument("-c", "--config", help="Yaml config file name", default="main")
    

    args = parser.parse_args()
    

    print(f"Protocol: {args.protocol}")
    print("------------------------")

    return args


def save_files(synapse, args, datestring, params):
    params = [x for x in params if len(x) > 0]
    fstr = ".".join(params)
    folder = f"{args.protocol}/{datestring}/{fstr}"
    synapse.write_debug_files(folder)


def run_simulation(function, datestring, args, namesubstr = "", neuron_parameters={}, synapse_parameters={}, alzheimers_parameters=[0.0, 0.0], beta_amyloid_parameters=[0.0, 1.0], configname=None):
    if float(alzheimers_parameters[0]) > 0:
        hocfile = "cell_seed4_0_AICD_cell.hoc"
    else:
        hocfile = "cell_seed4_0_control_cell.hoc"
    neuron_model = NeuronBuilder(hocfile=hocfile, recompile=False)
    synapse = SynapseType()

    neuron_model.set_parameters(neuron_parameters)
    synapse.set_parameters(synapse_parameters)

    neuron_model.set_alzheimers(*alzheimers_parameters)
    neuron_model.set_beta_amyloid(*beta_amyloid_parameters)

    getattr(neuron_model, function)(synapse)
    neuron_model.run()

    synapse.print_epsp_change()

    save_files(synapse, args, datestring, params = [ namesubstr ])
    
    if configname:
        config = load_config(configname['load'])
        save_config(config, f"debug/configs/{configname['save']}.yaml")


def all(args, datestring):
    processes = []

    NeuronBuilder.compile_mod_files(os.getcwd())
    functions = [
        "build_protocol_ltp",
        "build_protocol_ltd",
    ]

    configname = args.config
    config = load_config(configname)
    runs = config['runs']

    for runname in runs:
        runs[runname]['neuron_parameters'] = runs[runname]['neuron_parameters'] if 'neuron_parameters' in runs[runname] else {}
        runs[runname]['synapse_parameters'] = runs[runname]['synapse_parameters'] if 'synapse_parameters' in runs[runname] else {}
        for function in functions:
            namesubstr = function.split('_')[-1]
            p = Process(target=run_simulation, args=(
                function,
                f"{datestring}.{runname}",
                args,
                namesubstr,
                runs[runname]['neuron_parameters'],
                runs[runname]['synapse_parameters'],
                runs[runname]['alzheimers'],
                runs[runname]['beta_amyloid'],
                {'load': configname, 'save': datestring},
                ))
            p.start()
            processes.append(p)

    for p in processes:
        p.join()

def main():
    if 'MILEDIDEBUG' not in os.environ:
        os.environ['MILEDIDEBUG'] = '0'
    

    args = parse_args()
    datestring = datetime.now().strftime("%Y.%m.%d.%H.%M.%S")
    globals()[args.protocol](args, datestring)

if __name__ == "__main__":
    main()
