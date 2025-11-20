#!/usr/bin/env bash
# project/filestream.sh
set -euxo pipefail

# Export gwlandscape container SAS token (compasdatasets storage account)
export AZ_SAS="$GWLANDSCAPE_AZURE_BLOB_SAS"

# Full path for files list
FILE_LIST="/workspaces/cosmic_integration_dasein/project/files.txt"

# Loop through each line in the list and stream the file into container
while read -r url blob_path; do
  echo "Uploading $url as $blob_path ..."
  curl -L "$url" \
    | azcopy copy "https://compasdatasets.blob.core.windows.net/gwlandscape/$blob_path?$AZ_SAS" \
        --from-to=PipeBlob \
        --check-length=true
done < "$FILE_LIST"
