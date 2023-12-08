# IntegraTRN-Docker
Docker Image Platform for the IntegraTRN package

## IntegraTRN

The R package [IntegraTRN](https://github.com/j-y26/IntegraTRN) is hosted on https://github.com/j-y26/IntegraTRN.

See the package GitHub page for more details.

## Setting up a container of the Docker image:

### Install Docker

Ensure the operation system you use supports virtual environment.

To install Docker, go to the [Docker official page](https://www.docker.com/products/docker-desktop/) to install Docker Desktop. 

### Setting up the container

```bash
# Navigate to the directory where you would like to store your project

# Pull and run the docker image
docker run -e PASSWORD=changeit -v ${pwd}:/home/rstudio/projects -p 8787:8787 jyang26/integra_trn:v0.1.0
```

In your browser, navigate to http://localhost:8787 and log in with username `rstudio` and the password you specified.
