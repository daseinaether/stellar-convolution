# .devcontainer/post-create.sh
set -euxo pipefail


# --- AzCopy (Ubuntu 24.04) --- (https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-install-linux-package?tabs=apt)
curl -sSL -O https://packages.microsoft.com/config/ubuntu/24.04/packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
rm packages-microsoft-prod.deb
sudo apt-get update
sudo apt-get install -y azcopy


# --- micromamba init ---
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$(micromamba shell hook -s bash)"

# Future shell persistance
{
    echo 'export MAMBA_ROOT_PREFIX="$HOME/micromamba"'
    echo 'eval "$(micromamba shell hook -s bash)"'
} >> "$HOME/.bashrc"

# Strict channel enforcement (conda-forge)
micromamba config set channel_priority strict

# Create/update env from environment.yml
ENV_NAME="sspc"
if [ -f environment.yml ]; then
    micromamba create -y -n "$ENV_NAME" -f environment.yml || micromamba install -y -n "$ENV_NAME" -f environment.yml
else
    micromamba create -y -n "$ENV_NAME" python=3.13 
fi

# Auto-activate micromamba in shell
echo "micromamba activate $ENV_NAME" >> "$HOME/.bashrc"
