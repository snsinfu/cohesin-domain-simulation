#!/bin/sh -eu

mkdir -p _outputs _logs

concurrency_limit="${1:-100}"
config_count="$(jq '.["meta"].config_count' < config_template.json)"
task_spec="0-$((config_count - 1))%${concurrency_limit}"

qsub -t "${task_spec}" run.sh
