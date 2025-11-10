# .devcontainer/post-create.sh
set -euxo pipefail

# ===== AzCopy Installation (Ubuntu 24.04) =====
# Useful for streaming web-hosted datasets directly into Azure Blob Storage Containers.
#
# Allows for greater project flexibility, by removing the need to always carry the 
# physical hard drive containing the ~2+ TB datasets, and having code directly interface 
# with the cloud-hosted data anywhere in the world on any device.

# Installation Docs: https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-install-linux-package?tabs=apt
curl -sSL -O https://packages.microsoft.com/config/ubuntu/24.04/packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
rm packages-microsoft-prod.deb
sudo apt-get update
sudo apt-get install -y azcopy


# ===== Micromamba Init =====
# Many population synthesis and post-processing code is written inside a conda-forge Python environment.
#
# Conda-forge is typically setup via the Miniforge3 installer; 
# Micromamba is chosen here instead due to its compatibility with container environments.

# Initialize micromamba in the current shell
export MAMBA_ROOT="$HOME/micromamba"
eval "$(micromamba shell hook -s bash)"

# Intialize micromamba in future shell recurrences
{
    echo 'export MAMBA_ROOT="$HOME/micromamba"'
    echo 'eval "$(micromamba shell hook -s bash)"'
} >> "$HOME/.bashrc"

# Enforce conda-forge channel
micromamba config set channel_priority strict

# Create/update sspc env from environment.yaml
ENV="sspc"
if [ -f environment.yaml ]; then
    micromamba create -y -n "$ENV" -f environment.yml || micromamba install -y -n "$ENV" -f environment.yml
else
    micromamba create -y -n "$ENV" python=3.13 
fi

# Auto-activate env in the current shell
echo "micromamba activate $ENV" >> "$HOME/.bashrc"
