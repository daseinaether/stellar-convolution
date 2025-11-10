# .devcontainer/post-create.sh
set -euxo pipefail

# ===== AzCopy Installation =====
# Useful for streaming web-hosted datasets directly into Azure Blob Storage Containers.
#
# Allows for greater project flexibility, by removing the need to always carry the 
# physical hard drive containing the ~2+ TB datasets, and having code directly interface 
# with the cloud-hosted data anywhere in the world on any device.
sudo apt-get update
sudo apt-get install -y azcopy
azcopy --version || true


# ===== Micromamba Init =====
# Many population synthesis and post-processing code is written inside a conda-forge Python environment.
#
# Conda-forge is typically setup via the Miniforge3 installer; 
# Micromamba is chosen here instead due to its compatibility with container environments.
export MAMBA_ROOT="$HOME/micromamba"
eval "$(micromamba shell hook -s bash)"
{
  echo 'export MAMBA_ROOT="$HOME/micromamba"'
  echo 'eval "$(micromamba shell hook -s bash)"'
} >> "$HOME/.bashrc"

micromamba config set channel_priority strict

# Parse for the env config file, accepting either .yaml or .yml 
for F in environment.yaml environment.yml; do
  if [ -f "$F" ]; then ENV_FILE="$F"; break; fi
done

# Create/update the sspc env, or default to base env if ENV_FILE not found.
ENV=sspc
if [ -n "${ENV_FILE:-}" ]; then
  micromamba create -y -n "$ENV" -f "$ENV_FILE" || \
  micromamba env update -n "$ENV" -f "$ENV_FILE"
else
  micromamba create -y -n "base" python=3.12
fi

# Auto-activate env on future shell recurrences
echo "micromamba activate $ENV" >> "$HOME/.bashrc"

# Install pre-commit hooks if pre-commit package was successfully installed into env
if command -v pre-commit >/dev/null 2>&1; then
  pre-commit install || true
fi
