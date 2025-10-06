# AWS Batch deployment example

This directory provides a minimal example of how to run the tornado ground scoring
pipeline on [AWS Batch](https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html).
It uses a container image that downloads a configuration file from Amazon S3,
executes `python -m src.run_pipeline`, and then pushes the rendered outputs back to S3.

## Architecture overview

1. **Amazon S3** stores configuration files, AOI GeoJSON, and the pipeline outputs.
2. **Amazon Elastic Container Registry (ECR)** hosts the container image built from
   `Dockerfile` in this directory.
3. **AWS Batch** runs the container image on a compute environment (Fargate or EC2).
4. **Amazon EventBridge** (optional) triggers the Batch job on storm days, for example
   by reacting to a daily schedule or NOAA event feed.

## Build and publish the container image

```bash
# From the repository root
aws ecr create-repository --repository-name tornado-ground-scoring
ECR_URI=$(aws ecr describe-repositories \
  --repository-names tornado-ground-scoring \
  --query 'repositories[0].repositoryUri' --output text)

aws ecr get-login-password | docker login --username AWS --password-stdin "$ECR_URI"

docker build -f infra/aws-batch-example/Dockerfile -t tornado-ground-scoring:latest .
docker tag tornado-ground-scoring:latest "$ECR_URI:latest"
docker push "$ECR_URI:latest"
```

## Prepare S3 inputs

Upload your configuration YAML (and any referenced AOI or storm catalog files) to S3:

```bash
aws s3 cp data/config.yaml s3://my-bucket/configs/tornado-config.yaml
aws s3 cp data/storm_reports_example.csv s3://my-bucket/configs/storm_reports_example.csv
aws s3 cp data/example_aoi.geojson s3://my-bucket/configs/example_aoi.geojson
```

Update the S3 URIs inside the uploaded `config.yaml` so that file paths point to the
objects you just uploaded.

## Create the Batch job definition

Edit `job-definition.json` and replace the account IDs, role ARNs, and S3 bucket names
with your own. Then register the job definition:

```bash
aws batch register-job-definition \
  --cli-input-json file://infra/aws-batch-example/job-definition.json
```

You also need a compute environment and job queue. A minimal Fargate setup can be
created with:

```bash
aws batch create-compute-environment \
  --compute-environment-name tornado-ce \
  --type MANAGED \
  --state ENABLED \
  --compute-resources '{
    "type": "FARGATE",
    "maxvCpus": 32,
    "subnets": ["subnet-abc123"],
    "securityGroupIds": ["sg-abc123"]
  }'

aws batch create-job-queue \
  --job-queue-name tornado-queue \
  --priority 1 \
  --compute-environment-order '[{"order":1,"computeEnvironment":"tornado-ce"}]'
```

Finally, submit a test job:

```bash
aws batch submit-job \
  --job-name tornado-ground-scoring-test \
  --job-queue tornado-queue \
  --job-definition tornado-ground-scoring
```

The container downloads the config from `CONFIG_S3_URI`, runs the pipeline, and
synchronizes `docs/` and `artifacts/` back to `EXPORT_S3_PREFIX`.

## Automate with EventBridge (optional)

Create a rule that triggers only on storm days. For example, use a daily cron rule
that invokes a Lambda function. The Lambda can inspect the National Weather Service
storm reports API, and when a qualifying storm is detected, submit the Batch job via
`aws batch submit-job`.

```python
import json
import os

import boto3
import requests

BATCH_QUEUE = os.environ["BATCH_QUEUE"]
BATCH_JOB_DEFINITION = os.environ["BATCH_JOB_DEFINITION"]
STORM_API = "https://api.weather.gov/alerts/active?event=Tornado%20Warning"

batch = boto3.client("batch")


def handler(event, context):
    response = requests.get(STORM_API, timeout=10)
    response.raise_for_status()
    alerts = response.json().get("features", [])
    if not alerts:
        return {"triggered": False}

    batch.submit_job(
        jobName="tornado-ground-scoring-auto",
        jobQueue=BATCH_QUEUE,
        jobDefinition=BATCH_JOB_DEFINITION,
    )
    return {"triggered": True, "alerts": len(alerts)}
```

Use the Lambda to set `CONFIG_S3_URI` dynamically if you maintain per-day
configuration files.

## IAM considerations

Grant the Batch job role permission to:

- `s3:GetObject` for the configuration bucket/prefix.
- `s3:PutObject` and `s3:DeleteObject` for the output prefix.
- `logs:CreateLogStream` and `logs:PutLogEvents` for Amazon CloudWatch Logs.

The execution role should allow pulling the ECR image and writing logs.

## Cost controls

- Use Spot Fargate or EC2 Spot instances to reduce compute cost.
- Restrict EventBridge rules so that only storm-qualifying days trigger jobs.
- Lifecycle policies on the S3 output bucket to expire stale artifacts.

## Environment variables

The container entrypoint recognises the following environment variables:

| Name | Required | Description |
| ---- | -------- | ----------- |
| `CONFIG_S3_URI` | Yes | `s3://` URI of the main pipeline config file. |
| `EXTRA_CONFIG_S3_URI` | No | `s3://` prefix that will be recursively copied into the container before execution (for AOIs, storm catalogs, etc.). |
| `EXPORT_S3_PREFIX` | No | `s3://` prefix to which `docs/` and `artifacts/` are synced after the run. |
| `EXPORT_CSV_PATH` | No | Local path (inside the container) where the optional CSV export should be written; combine with `EXPORT_S3_PREFIX` to sync it out. |

Adjust the job definition to set these variables for each run.
