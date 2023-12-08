
## Python and NEURON codes for the NMDA-based synaptic plasticity model

This model analyses altered hippocampal synaptic plasticity and its rescue under the Alzheimer's disease (AD) conditions, when the concentrations of AD-related peptides, such as the amyloid precursor protein intracellular domain (AICD) and amyloid beta (Aβ), are increased. The phenomenological NMDA receptor-based voltage-dependent model is used to model synaptic modifications at the CA3-CA1 synapses onto the multicompartmental CA1 pyramidal neuron. The modeling results show that partial blockade of Glu2NB-NMDAR-gated channel restores intrinsic excitability of a CA1 pyramidal neuron and rescues long-term potentiation in AICD and Aβ conditions. The model is implemented in Python and NEURON.

### Reference
Justinas J Dainauskas, Paola Vitale, Sebastien Moreno, Helene Marie, Michele Migliore and Ausra Saudargiene. Altered synaptic plasticity at hippocampal CA1–CA3 synapses in Alzheimer's disease: integration of amyloid precursor protein intracellular domain and amyloid beta effects into computational models. Frontiers in Computational Neuroscience 2023 doi.org/10.3389/fncom.2023.1305169

The code reproduces Fig2-Fig9. 


## Run the code

```bash
python -m "src.manual_run" [options]
```

```
usage: manual_run.py [-h] [-p {ltp,ltd}] [-aC AICDCHANNELS] [-aN AICDNR2B]
                     [-n2b NR2BMULT] [-b BETA] [-bM BETA_MIDPOINT] [-s SESSION]

Build and run neuron simulation

options:
  -h, --help            show this help message and exit
  -p {ltp,ltd}, --protocol {ltp,ltd}
                        Provide protocol name.
                         ltp - 1 burst 100Hz 1s 
                         ltd - 1Hz 500s
  -aC AICDCHANNELS, --aicdchannels AICDCHANNELS
                        Fraction of AICD in channels [0-1]
  -aN AICDNR2B, --aicdnr2b AICDNR2B
                        Fraction of AICD for GluNR2B-NMDAR [0-4]
  -n2b NR2BMULT, --nr2bmult NR2BMULT
                        Fraction of GluNR2B-NMDAR  
  -b BETA, --beta BETA  Fraction of Beta Amyloids [0-1]
  -bM BETA_MIDPOINT, --beta_midpoint BETA_MIDPOINT
                        Coefficient for midpoint of LTP function to model beta amyloid effect
  -s SESSION, --session SESSION
                        Session file

```

For automatic calculations set parameters in "configs" directory yaml file and run with
```bash
python -m "src.auto_run" [options]
```

```
usage: auto_run.py [-h] [-p {all}] [-c CONFIG]

Build and run neuron simulation

options:
  -h, --help            show this help message and exit
  -p {all}, --protocol {all}
                        Protocol to run
  -c CONFIG, --config CONFIG
                        Yaml config file name
```

To plot the results calculated by auto_run calculations run:
```bash
python -m "src.plot" [options]
```

```
usage: plot.py [-h] [-p {all}] [-i ID]

Plot runs

options:
  -h, --help            show this help message and exit
  -p {all}, --protocol {all}
                        What to plot
  -i ID, --id ID        Which experiment to run from latest
```

## Reproduce figures

```bash
python -m src.auto_run -p all -c main && python -m src.plot -p all -i 0
```