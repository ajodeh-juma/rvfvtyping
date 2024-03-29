FROM nfcore/base:dev
LABEL authors="@ajodeh-juma" \
      description="Docker image containing all software requirements for the rvfvtyping pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/rvfvtyping-env/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name rvfvtyping-env > rvfvtyping-env.yml
