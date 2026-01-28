#!/bin/sh -eu

job_name="$(basename "$(pwd)")"
instance_id="$1"
log_prefix="_logs/${job_name}-${instance_id}"

../run.sh "${instance_id}" 1> "${log_prefix}.o" 2> "${log_prefix}.e"
