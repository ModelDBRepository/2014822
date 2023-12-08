#!/usr/bin/python
from .neuron_builder import NeuronBuilder
from .synapse_builder_ampanmdasynplast import AMPANMDASynPlastBuilder
from neuron import gui
import os
import argparse
import shutil
from datetime import datetime

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
                            choices=['ltp', 'ltd'],
                            default="ltp",
                            help="R|Provide protocol name.\n"
                            " ltp - 1 burst 100Hz 1s\n"
                            " ltd - 1Hz 500s\n"
                            )
    parser.add_argument("-aC", "--aicdchannels", help="Fraction of AICD in channels [0-1]", default="0")
    parser.add_argument("-aN", "--aicdnr2b", help="Fraction of AICD for GluNR2B-NMDAR [0-4]", default="1")
    parser.add_argument("-n2b", "--nr2bmult", help="Fraction of GluNR2B-NMDAR  ", default="1")    
    parser.add_argument("-b", "--beta", help="Fraction of Beta Amyloids [0-1]", default="0")
    parser.add_argument("-bM", "--beta_midpoint", help="Coefficient for midpoint of LTP function to model beta amyloid effect", default="1")
    parser.add_argument("-s", "--session", help="Session file", default="main.ses")

    args = parser.parse_args()
    
    if args.session != "none" and not os.path.exists(f"sessions/{args.session}"):
        print(f"No session file: sessions/{args.session}. Exiting...")
        exit(1)

    print(f"Protocol: {args.protocol}")
    print(f"AICD Channels: {args.aicdchannels}")
    print(f"AICD GluN2B: {args.aicdnr2b}")
    print("------------------------")

    return args
    

def ltp(args, datestring):
    if float(args.aicdchannels) > 0:
        hocfile = "cell_seed4_0_AICD_cell.hoc"
    else:
        hocfile = "cell_seed4_0_control_cell.hoc"
    neuron_model = NeuronBuilder(hocfile=hocfile)
    synapse = SynapseType()

    neuron_parameters = { 'nr2bmult': float(args.nr2bmult)}
    synapse_parameters = {}
    neuron_model.set_parameters(neuron_parameters)
    synapse.set_parameters(synapse_parameters)

    neuron_model.set_alzheimers(float(args.aicdchannels), float(args.aicdnr2b))
    neuron_model.set_beta_amyloid(float(args.beta), float(args.beta_midpoint))

    neuron_model.build_protocol_ltp(synapse)
    
    if args.session != "none":
        print(f"Loading session: {args.session}")
        neuron_model.load_session(f"sessions/{args.session}")    

    neuron_model.run()
    synapse.print_epsp_change()

def ltd(args, datestring):
    if float(args.aicdchannels) > 0:
        hocfile = "cell_seed4_0_AICD_cell.hoc"
    else:
        hocfile = "cell_seed4_0_control_cell.hoc"
    neuron_model = NeuronBuilder(hocfile=hocfile)
    synapse = SynapseType()

    neuron_parameters = { 'nr2bmult': float(args.nr2bmult)}
    synapse_parameters = {}
    neuron_model.set_parameters(neuron_parameters)
    synapse.set_parameters(synapse_parameters)

    neuron_model.set_alzheimers(float(args.aicdchannels), float(args.aicdnr2b))
    neuron_model.set_beta_amyloid(float(args.beta), float(args.beta_midpoint))

    neuron_model.build_protocol_ltd(synapse)

    if args.session != "none":
        print(f"Loading session: {args.session}")
        neuron_model.load_session(f"sessions/{args.session}")    

    neuron_model.run()
    synapse.print_epsp_change()

def main():
    global neuron_model, synapse

    if 'MILEDIDEBUG' not in os.environ:
        os.environ['MILEDIDEBUG'] = '0'

    if os.path.isdir('x86_64'):
        shutil.rmtree('x86_64')

    args = parse_args()


    datestring = datetime.now().strftime("%Y.%m.%d.%H.%M.%S")
    globals()[args.protocol](args, datestring)



if __name__ == "__main__":
    main()