from glob import glob
from pprint import pprint
import numpy as np
import os
import pickle


def read_protocol(experiment_id, protocol):
    print(f"Reading debug {protocol} experiment id: {experiment_id}...")
    files = list(filter(os.path.isdir, glob(f"debug/{protocol}/*")))
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    files = [os.path.basename(x) for x in files]
    pprint(files)

    experiment_path = files[experiment_id]
    files = list(filter(os.path.isdir, glob(f"debug/{protocol}/{experiment_path}/*")))
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    files = [os.path.basename(x) for x in files]
    
    all_data = {}
    for parameter_path in files:
        data = read_data(f"{protocol}/{experiment_path}/{parameter_path}")
        all_data[parameter_path] = data
        print(parameter_path)
        pprint(data['parameters'])
    
    return all_data


def read_protocol_bydatestring(datestring, protocol):
    print(f"Reading debug {protocol} datestring: {datestring}...")
    files = list(filter(os.path.isdir, glob(f"debug/{protocol}/{datestring}")))
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    files = [os.path.basename(x) for x in files]
    pprint(files)

    experiment_path = files[0]
    files = list(filter(os.path.isdir, glob(f"debug/{protocol}/{experiment_path}/*")))
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    files = [os.path.basename(x) for x in files]
    
    all_data = {}
    for parameter_path in files:
        data = read_data(f"{protocol}/{experiment_path}/{parameter_path}")
        all_data[parameter_path] = data
        print(parameter_path)
        if os.environ["MILEDIDEBUG"] == '1':
            pprint(data['parameters'])
    
    return all_data


def read_data(path):
    ec_synapse = np.load(f"debug/{path}/ec_synapse.npy", allow_pickle=True).tolist()
    synapse = np.load(f"debug/{path}/synapse.npy", allow_pickle=True).tolist()
    soma = np.load(f"debug/{path}/soma.npy", allow_pickle=True).tolist()
    t = np.load(f"debug/{path}/t.npy", allow_pickle=True).tolist()
    epsp = np.load(f"debug/{path}/epsp.npy", allow_pickle=True).tolist()
    parameters = np.load(f"debug/{path}/parameters.npy", allow_pickle=True).tolist()
    dendrite = np.load(f"debug/{path}/dendrite.npy", allow_pickle=True).tolist()
    if os.path.exists(f"debug/{path}/spikecount.npy"):
        spikes = np.load(f"debug/{path}/spikecount.npy", allow_pickle=True).tolist()
    else:
        spikes = []
    
    return {
        "t": t,
        "soma": soma,
        "synapse": synapse,
        "ec_synapse": ec_synapse,
        "epsp": epsp,
        "parameters": parameters,
        "dendrite": dendrite,
        "spikes": spikes
    }

