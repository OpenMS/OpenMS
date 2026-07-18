#!/usr/bin/env bash
# Dispatches a workflow_dispatch event via `gh workflow run`, retrying on
# transient failures (e.g. intermittent 5xx errors from the GitHub API).
set -uo pipefail

workflow="$1"
ref="$2"
max_attempts=3
delay=10

for attempt in $(seq 1 "$max_attempts"); do
  if gh workflow run "$workflow" --ref "$ref"; then
    exit 0
  fi
  if [ "$attempt" -lt "$max_attempts" ]; then
    echo "Attempt $attempt/$max_attempts to dispatch '$workflow' (ref: $ref) failed, retrying in ${delay}s..."
    sleep "$delay"
    delay=$((delay * 2))
  fi
done

echo "All $max_attempts attempts to dispatch '$workflow' (ref: $ref) failed."
exit 1
