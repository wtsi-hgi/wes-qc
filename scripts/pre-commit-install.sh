#!/bin/bash

conda activate hail
pip install pre-commit
pip install black mypy

pre-commit sample-config > .pre-commit-config.yaml

pre-commit install

pip install types-PyYAML
