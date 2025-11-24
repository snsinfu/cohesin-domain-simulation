#PBS -N d01
#PBS -o _logs
#PBS -e _logs
#PBS -l nodes=1:fast

cd "${PBS_O_WORKDIR:-.}"

forward_logs="n"

while getopts l opt; do
    case "${opt}" in
    l)  forward_logs="y" ;;
    *)  exit 1
    esac
done

shift $((OPTIND - 1))

task_id="${PBS_ARRAYID:-${1:-0}}"
config_id="${task_id}"

case "${forward_logs}" in
y)  exec >"_logs/${task_id}.log" 2>&1 ;;
*)
esac

seed="${task_id}"
config="_configs/config-${config_id}.json"
output="_outputs/output-${config_id}.h5"

echo c "${config}"
echo o "${output}"

sleep $((task_id % 10)).$((RANDOM % 10))

../../src/simulation/main -s "${seed}" -o "${output}" "${config}"
