#PBS -N d09
#PBS -o _logs
#PBS -e _logs
#PBS -l nodes=1:fast

cd "${PBS_O_WORKDIR:-.}"

task_id="${PBS_ARRAYID:-${1:-0}}"
config_id="${task_id}"

seed="${task_id}"
config="_configs/config-${config_id}.json"
output="_outputs/output-${config_id}.h5"

sleep $((task_id % 10)).$((RANDOM % 10))

../../../src/simulation/main -s "${seed}" -o "${output}" "${config}"
