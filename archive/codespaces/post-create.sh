#!/usr/bin/env bash
# /workspaces/cosmic_integration_dasein/.devcontainer/post-create.sh
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


# ===== Create/update sspc environment in micromamba =====

# Initialize vars for env setup
ENV="sspc"
CONFIG="/workspaces/cosmic_integration_dasein/.devcontainer/environment.yaml"

# Create/update sspc environment
# if [ -n "$CONFIG" ]; then
#   # Create env from file, if missing; Otherwise, update in place.
#   micromamba create -y -n "$ENV" -f "$CONFIG"
# else
#   # Announce error if config file missing; exit env setup process.
#   echo 'The "environment.yaml" configuration file was not found!'
# fi

# ========================================================


# ===== sspc env post-setup =====

# Create the sspc environment in micromamba
micromamba create -y -n "$ENV" -f "$CONFIG"

# Auto-activate env for future shell recurrences (idempotent append)
echo "micromamba activate $ENV" >> "$HOME/.bashrc"

# Install and run pre-commit git hooks
micromamba run -n "$ENV" pre-commit install || true

# Register Jupyter kernel
micromamba run -n "$ENV" python -m ipykernel install --user \
  --name "$ENV" --display-name "Python ($ENV)" || true

# Clean the container image
micromamba clean -a -y || true

# Announce completion
echo "Post-creation process complete!"

# ===============================
