#!/usr/bin/env bash
# .devcontainer/post-create.sh
set -euxo pipefail

# ===== Setup SSPC Environment in Micromamba =====

# Initialize vars for sspc env setup
ENV="sspc"
CONFIG=".devcontainer/environment.yaml"

# Create the sspc env in micromamba
micromamba create -y -n "$ENV" -f "$CONFIG"

# Auto-activate sspc env for future shell recurrences
echo "micromamba activate $ENV" >> "$HOME/.bashrc"

# Install pre-commit git hooks
micromamba run -n "$ENV" pre-commit install || true

# Register Jupyter kernel
micromamba run -n "$ENV" python -m ipykernel install --user \
  --name "$ENV" --display-name "Python ($ENV)" || true

# ===============================
