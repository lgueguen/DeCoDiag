# DecoDiag: Compute diagnostics of gene families from the incoherences
  in ancestral genomes

From repository src, docker container is built with command:

```{sh}
docker build -t anal_dock . -f Docker/Dockerfile 
```

Get the docker:
```{sh}
docker pull lgueguen/DecoDiag
```

Docker is run through script (logged in as root):

```{sh}
Docker/launch_docker.bash [working directory [parameter file]]
```

If parameter file argument is void, window interface is opened to ask.

if working directory argument is void, actual directory is used as
working directory.




