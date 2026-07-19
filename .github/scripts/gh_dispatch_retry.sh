#!/usr/bin/env bash
# Dispatch a workflow_dispatch event via `gh workflow run`, retrying on
# transient GitHub API errors (e.g. intermittent HTTP 5xx responses).
set -u

workflow="$1"
ref="$2"
max_attempts=5

for attempt in $(seq 1 "$max_attempts"); do
  if gh workflow run "$workflow" --ref "$ref"; then
    exit 0
  fi
  if [ "$attempt" -eq "$max_attempts" ]; then
    echo "::error::Failed to dispatch workflow '$workflow' on ref '$ref' after $max_attempts attempts"
    exit 1
  fi
  delay=$((5 * attempt))
  echo "Attempt $attempt/$max_attempts to dispatch '$workflow' on ref '$ref' failed, retrying in ${delay}s..."
  sleep "$delay"
done
