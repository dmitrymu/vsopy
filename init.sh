#!/usr/bin/env bash
python3 -m venv .venv
source ./.venv/bin/activate
pip install -r requirements.txt
git config --local core.hooksPath .githooks/
