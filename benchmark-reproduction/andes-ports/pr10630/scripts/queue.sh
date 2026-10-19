# usage: queue.sh <jobfile>; each line: build arm dataset [flags]; stops on first failure
while read -r line; do
  [ -z "$line" ] && continue
  echo "START $line $(date +%T)"
  python3 /home/user/bench/scripts/run_one.py $line > /dev/null 2>> /home/user/bench/queue_err.log || { echo "FAIL $line"; exit 1; }
  echo "DONE $line $(date +%T)"
done < "$1"
echo QUEUE_DONE
