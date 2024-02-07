# OpenMS Project Docker Images

The following images are made available [as packages in this repository](https://github.com/orgs/OpenMS/packages).

| Image Name | Description |
|------------|-------------|
| `openms-library` | Contains the OpenMS library and runtime dependencies |
| `openms-tools` (alias: `openms-executables`) | Contains the OpenMS library, tools, and runtime dependencies |
| `openms-tools-thirdparty` | Contains the OpenMS library, tools, thirdparty tools and runtime dependencies for all |

Within each image, the library and tools can all be found in `/opt/OpenMS`.

## Tags

Images are tagged with the release version, for example: `3.1.0`. The `latest` tag corresponds to the images built from the `master` and `nightly` branches.

## Use

To pull images, you can use the following command, substituting the image name and tag as described above.

```shell
docker pull ghcr.io/openms/openms-tools:3.1.0
```

Here's an example on running an tool:

```shell
docker run -t --rm ghcr.io/openms/openms-tools:3.1.0 IsobaricAnalyzer -h
```

The above command should output the help-text for `IsobaricAnalyzer`, and the container will be removed after. If you want to work on files in the directory you're issuing the command from, you can mount the directory as volume like so:

```shell
# downloading a file so this example works
wget https://github.com/OpenMS/OpenMS/raw/develop/share/OpenMS/examples/BSA/BSA1.mzML
docker run -t --rm -v "$PWD:/data" ghcr.io/openms/openms-tools:3.1.0 FileInfo -in /data/BSA1.mzML
```

For ease of use, you can even alias the above command. On linux you can add the following to your  `~/.bash_aliases` or `.bashrc` files:

```shell
alias openms='docker run -t --rm -v "$PWD:/data" ghcr.io/openms/openms-tools:3.1.0'
```

which will make the command significantly shorter: `openms FileInfo -in /data/BSA1.mzML`

## Development Notes

### Optimization Tips

To explore what is happening in each layer of the docker image, [`dive`](https://github.com/wagoodman/dive) can come in very handy. For ease of use, you can use this alias: `alias dive='docker run -ti --rm  -v /var/run/docker.sock:/var/run/docker.sock wagoodman/dive'`, and run `dive ghcr.io/openms/openms-tools:latest`.

To get the total size of an image: `docker image inspect -f "{{ .Size }}" ghcr.io/openms/openms-tools:latest | numfmt --to=si`
