# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


# build bioconductor_rnaseq docker image

This assignment walks you through modifying an [bioconductor docker images](https://hub.docker.com/r/bioconductor/bioconductor_docker) by installing a variety of bioconductor packages as well as [asciinema](https://asciinema.org/docs/installation) and saving the new docker image to dockerhub. [asciinema](https://asciinema.org/about) is a free and open source solution for recording terminal sessions and sharing them on the web. [Docker](https://www.docker.com/resources/what-container) is tool that packages software into self-contained computing environments, called containers. [Docker Hub](https://www.docker.com/products/docker-hub) is a hosted repository service provided by Docker for finding and sharing container images with others.

<!-- blank line -->
----
<!-- blank line -->

## Learning Objectives:
 - pull a [bioconductor docker image](https://hub.docker.com/r/bioconductor/bioconductor_docker) from DockerHub
 - run the bioconductor/bioconductor_docker container via docker
 - install asciinema via CLI
 - save changes to a new docker image
 - push the new docker image to dockerhub
 
  ## Assignment 
1. Complete the assignment described below.
2. Upload a link to your dockerhub account.
3. Upload a link with screen-cast.

### Prerequisites
* create an asciinema account using email at [asccinema.org](https://asciinema.org/login/new) 
* create a [dockerhub account](https://hub.docker.com/)
* navigate to the directory: ~/assignments/bioconductor_asciinema
<!-- blank line -->
----
<!-- blank line -->

 ### Assignment Points
|  Rubric        | Points | 
|----------------|-------|
| Dockerhub     |  -/5  |
| Screencast     |  -/5  |
| On Time        |  -/5  |
*Total Points: -/10*

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

### 3. boot into the new docker container 
```
docker run -it bioconductor_asciinema:latest bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. start a screen-cast from within the container 
```
asciinema rec
```
<!-- blank line -->
----
<!-- blank line -->

### 5. link your container to your asciinema.org account by opening the URL in a web browser 
```
asciinema auth
```
<!-- blank line -->
----
<!-- blank line -->

### 6. add screen-cast headers and check your docker container
```
# GMS6804
# name: [your-name-here]
# date: [current-date]
# semester: [current-semester]
# assignment: build bioconductor-asciinema  

# check if R is installed?
R
# exit R
quit()
```
<!-- blank line -->
----
<!-- blank line -->

### 7. login to dockerhub
```
docker login
login: [YOUR DOCKERHUB ID]
pwd: [YOUR PASSWORD]
```
<!-- blank line -->
----
<!-- blank line -->

### 8. tag container
```
docker tag bioconductor_asciinema [YOUR DOCKERHUB ID]/bioconductor_asciinema:[month_year]
docker push [YOUR DOCKERHUB ID]/bioconductor_asciinema:[month_year]
```

### 9. example screen-cast