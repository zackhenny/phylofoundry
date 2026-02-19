#!/bin/bash
# Wrapper to run phylofoundry from source without installation
export PYTHONPATH=$PYTHONPATH:$(pwd)/src
python3 -m phylofoundry.main "$@"
