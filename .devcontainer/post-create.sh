# .devcontainer/post-create.sh
set -euxo pipefail


# Install AzCopy (https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-install-linux-package?tabs=apt)
curl -sSL -O https://packages.microsoft.com/config/ubuntu/24.04/packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
rm packages-microsoft-prod.deb
sudo apt-get update
sudo apt-get install -y azcopy


# Enable micromamba via shell hook.
echo 'eval "$(micromamba shell hook -s bash)"' >> "$HOME/.bashrc"
source "$HOME/.bashrc"

# Create sspc env from file. Defaults to base env creation if missing file.
if [ -f environment.yml ]; then
    micromamba create -y -n sspc -f environment.yml
    echo 'micromamba activate sspc' >> "$HOME/.bashrc"
else
    micromamba create -y -n base python=3.14
    echo 'micromamba activate base' >> "$HOME/.bashrc"
fi
