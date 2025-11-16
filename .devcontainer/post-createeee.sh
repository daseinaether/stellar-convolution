#!/usr/bin/env bash
# .devcontainer/post-create.sh
set -euxo pipefail

# ===== AzCopy Installation =====
# Useful for streaming web-hosted datasets directly into Azure Blob Storage Containers.
#
# Allows for greater project flexibility, by removing the need to always carry the
# physical hard drive containing the ~2+ TB datasets, and having code directly interface
# with the cloud-hosted data anywhere in the world on any device.

# Download repository configuration package
curl -sSL -O https://packages.microsoft.com/config/ubuntu/24.04/packages-microsoft-prod.deb

# Install repository configuration package
sudo dpkg -i packages-microsoft-prod.deb

# Remove repository configuration package after installation
rm packages-microsoft-prod.deb

# Update package index files
sudo apt-get update

# Install AzCopy
sudo apt-get install azcopy

# Announce AzCopy version if installed successfully
azcopy --version || true
# ===============================

# # # Micromamba is already initialized during the container build process
# # ===== Micromamba Init =====
# # Many stellar population synthesis and post-processing code
# # is written inside a conda-forge Python environment.
# #
# # Conda-forge is typically setup via the Miniforge3 installer;
# # Micromamba is chosen here instead due to its compatibility with container environments.

# # Setup shell hook for current instance
# export MAMBA_ROOT="$HOME/micromamba"
# eval "$(micromamba shell hook -s bash)"

# # Persist shell hook for future recurrences (idempotent append)
# grep -q 'micromamba shell hook -s bash' "$HOME/.bashrc" || {
#   echo 'export MAMBA_ROOT="$HOME/micromamba"' >> "$HOME/.bashrc"
#   echo 'eval "$(micromamba shell hook -s bash)"' >> "$HOME/.bashrc"
# }

# # The micromamba feature already defaults to strict channel_priority (conda-forge)
# micromamba config set channel_priority strict

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

# Install and run pre-commit git hooks
micromamba run -n "$ENV" pre-commit install >/dev/null 2>&1 || true
micromamba run -n "$ENV" pre-commit run --all-files >/dev/null 2>&1 || true

# Register Jupyter kernel
micromamba run -n "$ENV" python -m ipykernel install --user \
  --name "$ENV" --display-name "Python ($ENV)" >/dev/null 2>&1 || true

# Clean the container image
micromamba clean -a -y || true

# Announce completion
echo "Post-creation process complete!"

# ===============================
