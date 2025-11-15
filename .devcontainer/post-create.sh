# .devcontainer/post-create.sh
#!/usr/bin/env bash
set -euxo pipefail

# ===== AzCopy Installation =====
# Useful for streaming web-hosted datasets directly into Azure Blob Storage Containers.
#
# Allows for greater project flexibility, by removing the need to always carry the
# physical hard drive containing the ~2+ TB datasets, and having code directly interface
# with the cloud-hosted data anywhere in the world on any device.

export DEBIAN_FRONTEND=noninteractive
sudo apt-get update
sudo apt-get install -y azcopy
azcopy --version || true
# ===============================


# ===== Micromamba Init =====
# Many stellar population synthesis and post-processing code
# is written inside a conda-forge Python environment.
#
# Conda-forge is typically setup via the Miniforge3 installer;
# Micromamba is chosen here instead due to its compatibility with container environments.

# Setup shell hook for current instance
export MAMBA_ROOT="$HOME/micromamba"
eval "$(micromamba shell hook -s bash)"

# Persist shell hook for future recurrences (idempotent append)
grep -q 'micromamba shell hook -s bash' "$HOME/.bashrc" || {
  echo 'export MAMBA_ROOT="$HOME/micromamba"' >> "$HOME/.bashrc"
  echo 'eval "$(micromamba shell hook -s bash)"' >> "$HOME/.bashrc"
}

# Deterministic channel priority resolution
micromamba config set channel_priority strict

# Initialize env setup vars
ENV="sspc"
ENV_FILE=""

# Parse for env config file in root dir
for F in environment.yml environment.yaml; do
  if [ -f "$F" ]; then
    ENV_FILE="$F"
    break
  fi
done

# Create/update sspc environment
if [ -n "$ENV_FILE" ]; then
  # Create from file if missing; Otherwise, update in place.
  micromamba create -y -n "$ENV" -f "$ENV_FILE" || \
  micromamba env update -n "$ENV" -f "$ENV_FILE"
else
  # Fallback: Default to Python 3.12 (same as pyproject.toml)
  micromamba create -y -n "$ENV" python=3.12
fi
# ===========================


# ===== sspc env post-setup =====

# Auto-activate env for future shell recurrences (idempotent append)
grep -q "micromamba activate $ENV" "$HOME/.bashrc" || \
echo "micromamba activate $ENV" >> "$HOME/.bashrc"

# Install pre-commit git hooks (if pre-commit package successfully installed)
if micromamba run -n "$ENV" pre-commit --version >/dev/null 2>&1; then
  micromamba run -n "$ENV" pre-commit install || true
fi

# Register Jupyter kernel (if ipykernel package successfully installed)
micromamba run -n "$ENV" python -m ipykernel install --user --name "$ENV" --display-name "Python ($ENV)" || true

# Clean the container image
micromamba clean -a -y || true

# Announce completion
echo "Post-creation process complete!"

# ===============================
