FROM mambaorg/micromamba:1.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml

# Install dependencies
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Set working directory
WORKDIR /app

# Copy package source
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Install the package
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install .

# Ensure correct permissions
USER $MAMBA_USER

# Set entrypoint to the installed CLI tool
ENTRYPOINT ["phylofoundry"]
CMD ["--help"]
