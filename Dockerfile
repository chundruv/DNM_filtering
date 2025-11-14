# Multi-stage build for smaller, secure images
# Build stage - contains compilation tools
FROM python:3.11-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy and install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Runtime stage - minimal dependencies
FROM python:3.11-slim

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Create non-root user
RUN useradd -m -s /bin/bash optimizer
USER optimizer
WORKDIR /home/optimizer

# Copy package
COPY --chown=optimizer:optimizer src/ ./src/
COPY --chown=optimizer:optimizer pyproject.toml .

# Install package in development mode
RUN pip install -e .

# Create directories for data and results
RUN mkdir -p data results cache

# Default command
CMD ["variant-optimize", "--help"]

# Usage:
# docker build -t variant-optimizer .
# docker run -v $(pwd)/data:/home/optimizer/data \
#            -v $(pwd)/results:/home/optimizer/results \
#            variant-optimizer \
#            variant-optimize data/variants.tsv data/reference.tsv \
#            --output results --preset balanced
