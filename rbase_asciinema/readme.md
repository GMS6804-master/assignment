# build asciinema-rbase docker image

This assignment walks you through modifying an [r-base docker image](https://hub.docker.com/_/r-base?tab=description) by installing [asciinema](https://asciinema.org/docs/installation) and saving the new docker image to dockerhub.[asciinema](https://asciinema.org/about) is a free and open source solution for recording terminal sessions and sharing them on the web. [Docker](https://www.docker.com/resources/what-container) is tool that packages software into self-contained computing environments, called containers. [Docker Hub](https://www.docker.com/products/docker-hub) is a hosted repository service provided by Docker for finding and sharing container images with others.
<!-- blank line -->
----
<!-- blank line -->
## Learning Objectives:
 - pull a [r-base docker image](https://hub.docker.com/_/r-base?tab=description) from DockerHub
 - download ~/assignment/repository from Github
 - run the r-base container via docker
 - create a new container that modifies the r-base image using a dockerfile 
 - save new docker container to dockerhub
 - upload your screen-cast recording(s) to asciinema.org
 
  ## Assignment 
1. Complete the assignment described below.
2. Upload a link with screen-cast to the CANVAS Assignment.

### Prerequisites
* create an asciinema account using email at [asccinema.org](https://asciinema.org/login/new) 
* create a [dockerhub account](https://hub.docker.com/)
* download the ~/assignments/ GitHub repository. 
* navigate to the directory: ~/assignments/build_rbase_asciinema
<!-- blank line -->
----
<!-- blank line -->

## Getting Started

### 1. open docker teminal
<!-- blank line -->
----
<!-- blank line -->

### 2. pull the docker image
```
docker pull r-base
```
<!-- blank line -->
----
<!-- blank line -->

### 3. boot into a docker container 
```
docker run -it r-base:latest bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. start a screen-cast recording- what is the problem? 
```
asciinema rec
```
<!-- blank line -->
----
<!-- blank line -->

### 5. create a new container that installs asciinema using a dockerfile
```
docker build -t rbase_asciinema .
```
<!-- blank line -->
----
<!-- blank line -->

### 6. boot into a docker container. this requires that your working directory is ~/assignments/build_rbase_asciinema
```
docker run -it rbase_asciinema:latest bash
```
<!-- blank line -->
----
<!-- blank line -->

### 7. start a screen-cast from within the container 
```
asciinema rec
```
<!-- blank line -->
----
<!-- blank line -->

### 8. link your container to your asciinema.org account by opening the URL in a web browser 
```
asciinema auth
```
<!-- blank line -->
----
<!-- blank line -->

### 9. add screen-cast headers and modify your docker container
```
# GMS6804
# name: [your-name-here]
# date: [current-date]
# semester: [current-semester]
# assignment: build rbase-asciinema  

# check if R is installed?
R
# exit R
quit()
# exit container 
ctr-d
```
<!-- blank line -->
----
<!-- blank line -->

### 10. login to dockerhub
```
docker login
login: [YOUR DOCKERHUB ID]
pwd: [YOUR PASSWORD]
```
<!-- blank line -->
----
<!-- blank line -->

### 11. tag container
```
docker tag rbase_asciinema [YOUR DOCKERHUB ID]/rbase_asciinema:[month_year]
docker push [YOUR DOCKERHUB ID]/rbase_asciinema:[month_year]
```

### 12. example screen-cast
