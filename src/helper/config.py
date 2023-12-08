import yaml
from pathlib import Path

def load_debug_config(configname):
    with open(f"debug/configs/{configname}.yaml", 'r') as stream:
        data_loaded = yaml.safe_load(stream)
    return data_loaded

def load_config(configname):
    with open(f"configs/{configname}.yaml", 'r') as stream:
        data_loaded = yaml.safe_load(stream)
    return data_loaded

def save_config(config_data, filepath):
    Path(filepath).parent.mkdir(exist_ok = True, parents = True)
    with open(filepath, 'w') as stream:
        yaml.safe_dump(config_data, stream)
