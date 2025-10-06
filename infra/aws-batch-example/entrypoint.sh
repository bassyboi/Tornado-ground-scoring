#!/bin/bash
set -euo pipefail

if [[ -z "${CONFIG_S3_URI:-}" ]]; then
  echo "CONFIG_S3_URI environment variable must be set to an s3:// URI" >&2
  exit 1
fi

WORKDIR="/workspace"
CONFIG_PATH="$WORKDIR/config.yaml"

mkdir -p "$WORKDIR/artifacts" "$WORKDIR/docs"

aws s3 cp "$CONFIG_S3_URI" "$CONFIG_PATH"

if [[ -n "${EXTRA_CONFIG_S3_URI:-}" ]]; then
  aws s3 cp "$EXTRA_CONFIG_S3_URI" "$WORKDIR/" --recursive
fi

python3 -m src.run_pipeline --config "$CONFIG_PATH" ${EXPORT_CSV_PATH:+--export-csv "$EXPORT_CSV_PATH"}

if [[ -n "${EXPORT_S3_PREFIX:-}" ]]; then
  aws s3 sync "$WORKDIR/docs" "${EXPORT_S3_PREFIX%/}/docs" --delete
  aws s3 sync "$WORKDIR/artifacts" "${EXPORT_S3_PREFIX%/}/artifacts" --delete
fi
