# build bioconductor_rnaseq docker image

This assignment walks you through modifying an [bioconductor docker images](https://hub.docker.com/r/bioconductor/bioconductor_docker) by installing a variety of bioconductor packages as well as [asciinema](https://asciinema.org/docs/installation) and saving the new docker image to dockerhub. [asciinema](https://asciinema.org/about) is a free and open source solution for recording terminal sessions and sharing them on the web. [Docker](https://www.docker.com/resources/what-container) is tool that packages software into self-contained computing environments, called containers. [Docker Hub](https://www.docker.com/products/docker-hub) is a hosted repository service provided by Docker for finding and sharing container images with others.

<!-- blank line -->
----
<!-- blank line -->

## Learning Objectives:
 - build a [bioconductor docker image](https://hub.docker.com/r/bioconductor/bioconductor_docker) from a dockerfile
 - push the new docker image to dockerhub
 
  ## Assignment 
1. Complete the assignment described below.
2. Upload the link to your dockerhub account.


### Prerequisites
* create a [dockerhub account](https://hub.docker.com/)
* navigate to the directory: ~/path/to/directory/with/dockerfile/
<!-- blank line -->
----
<!-- blank line -->

 ### Assignment Points
|  Rubric        | Points | 
|----------------|-------|
| Dockerhub     |  -/50  |
| On Time        |  -/20  |
*Total Points: -/70*

## Getting Started

### 1. open docker teminal
<!-- blank line -->
----
<!-- blank line -->

### 2. create a new container that installs asciinema using a dockerfile
```
docker build -t bioconductor_asciinema .
```
<!-- blank line -->
----
<!-- blank line -->

### 3. test your the new docker container 
```
docker run -it bioconductor_asciinema:latest bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. login to dockerhub
```
docker login
login: [YOUR DOCKERHUB ID]
pwd: [YOUR PASSWORD]
```
<!-- blank line -->
----
<!-- blank line -->

### 5. tag container
```
docker tag bioconductor_asciinema [YOUR DOCKERHUB ID]/bioconductor_asciinema:[month_year]
```
<!-- blank line -->
----
<!-- blank line -->

### 6. push to DockerHub
```
docker push [YOUR DOCKERHUB ID]/bioconductor_asciinema:[month_year]
```
