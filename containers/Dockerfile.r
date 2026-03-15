# containers/Dockerfile.r
FROM bioconductor/bioconductor_docker:RELEASE_3_18

WORKDIR /app

COPY envs/env_r.yaml /tmp/env_r.yaml

# Install conda for R env management
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda env create -f /tmp/env_r.yaml && \
    /opt/conda/bin/conda clean -afy

# Install GitHub-only R packages
RUN /opt/conda/envs/env_r/bin/Rscript \
    -e "remotes::install_github('nalab-stanford/ALRA', upgrade='never')" \
    -e "remotes::install_github('IanevskiAleksandr/sc-type', upgrade='never')"

COPY src/scanpy_workflow/ambient/   /app/ambient/
COPY src/scanpy_workflow/preprocessing/ /app/preprocessing/
COPY src/scanpy_workflow/imputation/    /app/imputation/
COPY src/scanpy_workflow/annotation/    /app/annotation/
