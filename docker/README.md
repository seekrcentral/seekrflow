You need to create a personal access token from Github and place it into the
directory: ~/.config/seekr/github_token.

Grant access to private seekr repository: 
* Contents - read only.

NOTE: all commands must be run from seekrflow/ base directory.

# Build (from the seekrflow repo root):

   DOCKER_BUILDKIT=1 docker build \
  --network=host \
  -f docker/seekr-base.Dockerfile \
  --secret id=github_token,src=$HOME/.config/seekr/github_token \
  --build-arg SEEKR_CACHE_BUST="$(date +%s)" \
  -t seekr-base:latest .

# Run the following command to see the outputs from the previous run:

docker run --rm seekr-base:latest \
    python -c "import openmm; print('OpenMM', openmm.version.version); print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"

# Next, test out the CUDA platform

docker run --rm --gpus all seekr-base:latest \
    python -m openmm.testInstallation

# Next, Build the BD engine Dockerfile

DOCKER_BUILDKIT=1 docker build \
     --network=host \
     -f docker/seekr-engines-bd.Dockerfile \
     -t seekr-engines-bd:latest .

# Run the following command to see the outputs from the previous run:

docker run --rm seekr-engines-bd:latest \
    which bd_top || echo "WARNING: bd_top not found on PATH"