#PBS -N BP3
#PBS -o _logs
#PBS -e _logs
#PBS -l nodes=1:fast

cd "${PBS_O_WORKDIR:-.}"

task_id="${PBS_ARRAYID:-${1:-0}}"

replica_count="${REPLICA_COUNT:-10}"
replica_offset="${REPLICA_OFFSET:-0}"
config_id="$((task_id / replica_count))"
instance_id="$((task_id % replica_count + replica_offset))"

seed_offset="${SEED_OFFSET:-0}"
seed="$((task_id + seed_offset))"
config="_configs/config-${config_id}.json"
output="_outputs/output-${config_id}-${instance_id}.h5"

# Desync I/O
sleep $((task_id % 10)).$((RANDOM % 10))

../../../src/simulation/main -s "${seed}" -o "${output}" "${config}"
