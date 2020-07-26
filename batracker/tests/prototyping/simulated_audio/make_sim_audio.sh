#!/bin/bash
conda init bash
conda activate batracker
python -m tacost -paramfile sim_audio_params.yaml
