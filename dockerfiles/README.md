# OpenMS Project Docker Images

The following images are made available:

| Image Name | Description |
|------------|-------------|
| `openms-library` | Contains the OpenMS library and runtime dependencies |
| `openms-tools` | Contains the OpenMS library, tools, and runtime dependencies |
| `openms-tools-thirdparty` | Contains the OpenMS library, tools, thirdparty tools and runtime dependencies for all |

Within each image, the library and tools can all be found in `/opt/OpenMS`.

## Development Notes

### Optimization Tips

To explore what is happening in each layer of the docker image, [`dive`](https://github.com/wagoodman/dive) can come in very handy. For ease of use, you can use this alias: `alias dive='docker run -ti --rm  -v /var/run/docker.sock:/var/run/docker.sock wagoodman/dive'`, and run `dive ghcr.io/openms/openms-tools:latest`.

To get the total size of an image: `docker image inspect -f "{{ .Size }}" ghcr.io/openms/openms-tools:latest | numfmt --to=si`
