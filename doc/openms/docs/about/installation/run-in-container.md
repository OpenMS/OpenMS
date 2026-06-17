# Container

% start-after
## Run OpenMS inside a (Bio)Container

1. Install a containerization software (e.g., [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps))

2. Pull an image from one of the following registries:

- [OpenMS GitHub Container Registry](https://ghcr.io) for nightly binaries AND releases:

	On our registry, we provide one image for the library (with contrib) and one for the executables (with thirdparty).

	1. [openms-library](https://ghcr.io/openms/openms-library)
	2. [openms-executables](https://ghcr.io/openms/openms-executables)

	They can be pulled/run via the following commands:

	::::{tab-set}

	:::{tab-item} Docker

	```{code-block} bash
	docker pull ghcr.io/openms/openms-library
	docker pull ghcr.io/openms/openms-executables
	```

	:::

	:::{tab-item} Singularity

	```{code-block} bash
	singularity run ghcr.io/openms/openms-library-sif
	singularity run ghcr.io/openms/openms-executables-sif
	```

	:::

	::::

	:::{note}
	Per default this results in the download of the latest nightly snapshot. Specify a release version (e.g.,
	{{ '`docker pull ghcr.io/openms/openms-library:{0}`'.format(version) }} to receive a stable version.
	:::

- Otherwise, the [BioContainers Registries](https://biocontainers.pro/registry) and the associated Galaxy
project provide native containers based on our bioconda packages for both Docker and Singularity.

	1. [BioContainers libopenms](https://biocontainers.pro/tools/libopenms)
	2. [BioContainers openms](https://biocontainers.pro/tools/openms)
	3. [BioContainers openms-thirdparty](https://biocontainers.pro/tools/openms-thirdparty)
	4. [BioContainers pyOpenMS](https://biocontainers.pro/tools/pyopenms)

	Images of the containers can be pulled via or one of the following commands:

	::::{tab-set}

	:::{tab-item} Docker

	```{code-block} bash
	docker pull quay.io/biocontainers/libopenms
	docker pull quay.io/biocontainers/openms
	docker pull quay.io/biocontainers/pyopenms
	docker pull quay.io/biocontainers/openms-thirdparty
	```

	:::

	:::{tab-item} Singularity

	```{code-block} bash
	singularity run https://depot.galaxyproject.org/singularity/libopenms
	singularity run https://depot.galaxyproject.org/singularity/openms
	singularity run https://depot.galaxyproject.org/singularity/pyopenms
	singularity run https://depot.galaxyproject.org/singularity/openms-thirdparty
	```

	:::

	::::

:::{note}
If Singularity images fail to download or run, try to use the Docker images as Singularity will automatically convert them.
:::

Dockerfiles to build different kind of images (e.g., for ArchLinux) yourself can be found on
GitHub in our [OpenMS/dockerfiles](https://github.com/OpenMS/dockerfiles) repository. They usually follow our build
instructions closely, so you can have a look on how this is done in a clean environment.