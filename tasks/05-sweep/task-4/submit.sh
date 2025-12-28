#!/bin/sh -eu

sweep_count="$(jq '.["@meta"].sweep_count' < config_template.json)"
max_concurrency=100
array_spec="0-$((sweep_count - 1))%${max_concurrency}"

qsub -N "d05.4" -t "${array_spec}" ../run.sh
