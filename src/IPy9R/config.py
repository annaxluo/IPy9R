# src/IPy9R/config.py

import yaml


def load_config(config_fn):
    with open(config_fn, "r") as f:
        return yaml.safe_load(f)


def get_step_config(config, step_name):
    shared = config.get("shared", {})
    step = config.get(step_name, {})
    return {**shared, **step}