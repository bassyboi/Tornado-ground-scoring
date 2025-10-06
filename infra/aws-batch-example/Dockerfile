# Example container for running the tornado ground-scoring pipeline on AWS Batch.
# Build with: docker build -f infra/aws-batch-example/Dockerfile .

FROM ghcr.io/osgeo/gdal:ubuntu-small-3.7.3

WORKDIR /workspace

# Install Python dependencies and AWS CLI.
COPY requirements.txt ./
RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install --no-cache-dir -r requirements.txt awscli

# Copy the project code into the image.
COPY src ./src
COPY data ./data
COPY docs ./docs
COPY infra/aws-batch-example/entrypoint.sh ./infra/aws-batch-example/entrypoint.sh
COPY README.md ./README.md

ENV PYTHONPATH=/workspace

ENTRYPOINT ["/workspace/infra/aws-batch-example/entrypoint.sh"]
