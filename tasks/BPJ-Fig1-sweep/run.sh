#PBS -N BP1s
#PBS -o _logs
#PBS -e _logs
#PBS -l nodes=1:fast

cd "${PBS_O_WORKDIR:-.}"

task_id="${PBS_ARRAYID:-${1:-0}}"

instance_count=1
config_id="$((task_id / instance_count))"
instance_id="$((task_id % instance_count))"

seed="${task_id}"
config="_configs/config-${config_id}.json"
output="_outputs/output-${config_id}-${instance_id}.h5"

# Desync I/O
sleep $((task_id % 10)).$((RANDOM % 10))

../../../src/simulation/main -s "${seed}" -o "${output}" "${config}"
