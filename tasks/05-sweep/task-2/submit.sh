#!/bin/sh

sweep_count="$(jq '.["@meta"].sweep_count' < config_template.json)"
max_concurrency=300

qsub -t "0-${sweep_count}%${max_concurrency}" ../run.sh
