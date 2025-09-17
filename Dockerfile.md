# Dockerfile for Glimepiride Model

## Build image
To build the latest image use:
```bash
docker build -f Dockerfile -t matthiaskoenig/glimepiride:1.1 -t matthiaskoenig/glimepiride:latest .
```

## Push images
The image is pushed to dockerhub: [Docker Hub – Glimepiride](https://hub.docker.com/repository/docker/matthiaskoenig/glimepiride/general)

```
docker login
docker push --all-tags matthiaskoenig/glimepiride
```

## Run container
To use the latest container version interactively use:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/glimepiride:latest /bin/bash
```

To use a specific container version provide the version tag:
```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/glimepiride:1.1 /bin/bash
```

## Run simulations
Run the complete analysis:
```bash
uv run run_glimepiride -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```
