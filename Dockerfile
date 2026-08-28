# REDQuanTEA Docker Image
# Replication-Enhanced Detection of Quantitative Traits Evolving Adaptively
#
# Build: docker build -t redquantea .
# Run Detection Module: docker run -v $(pwd)/results:/app/results redquantea snakemake --configfile config/config_detect.yaml --cores 4
# Run Design Module: docker run -v $(pwd)/results:/app/results redquantea snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 2

FROM mambaforge/mambaforge:latest

LABEL maintainer="REDQuanTEA Team"
LABEL description="REDQuanTEA: Replication-Enhanced Detection of Quantitative Traits Evolving Adaptively"
LABEL version="1.0.0"

# Set working directory
WORKDIR /app

# Copy environment file first (for better caching)
COPY environment.yml /app/

# Create conda environment
RUN mamba env create -f environment.yml && \
    mamba clean --all --yes

# Make conda environment default
SHELL ["conda", "run", "-n", "redquantea", "/bin/bash", "-c"]

# Copy repository contents
COPY . /app/

# Create results directory
RUN mkdir -p /app/results

# Set environment variables
ENV PATH="/opt/conda/envs/redquantea/bin:$PATH"
ENV CONDA_DEFAULT_ENV=redquantea

# Verify installation
RUN Rscript -e "library(abc); library(ggplot2); cat('R packages verified\n')" && \
    python -c "import snakemake; print('Snakemake verified')"

# Default command: show help
CMD ["snakemake", "--help"]

# Entry point for running Snakemake workflows
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "redquantea"]
