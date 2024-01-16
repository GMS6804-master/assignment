# my first screen-cast 

This assignment walks you through creating a screen-cast via asciinema and docker. [asciinema](https://asciinema.org/about) is a free and open source solution for recording terminal sessions and sharing them on the web. [Docker](https://www.docker.com/resources/what-container) is tool that packages software into self-contained computing environments, called containers.
<!-- blank line -->
----
<!-- blank line -->
## Learning Objectives:
 - pull a [asciinema docker image](https://hub.docker.com/r/asciinema/asciinema/) from DockerHub
 - run the asciinema/asciinema container via docker
 - create a screen-cast recording
 - upload your recording to asciinema.org
 
 ## Assignment 
1. Complete the assignment described below.
2. Upload a link with screen-cast to the CANVAS Assignment.
 
 ### Assignment Points
|  Rubric        | Points | 
|----------------|-------|
| Screencast     |  -/5  |
| On Time        |  -/5  |

*Total Points: -/10*

### Prerequisites
* create and account using email with [asccinema.org](https://asciinema.org/login/new) 
<!-- blank line -->
----
<!-- blank line -->

## Getting Started
### 1. open docker teminal

![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/terminal_start.png)
<!-- blank line -->
----
<!-- blank line -->

### 2. pull the docker image

```
docker pull dominicklemas/rbase_asciinema:12_2021
```
<!-- blank line -->
----
<!-- blank line -->

### 3.  create a docker container 

```
docker run -it dominicklemas/rbase_asciinema:12_2021 bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. start a screen recording session

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
![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/asciinema_auth.png)
<!-- blank line -->
----
<!-- blank line -->

### 6. create your first screen-cast recording
```
# hello translational bioinformatics!
# my name is [your-name-here]
# todays date is [current-date]

```
<!-- blank line -->
----
<!-- blank line -->

### 7. stop your screen-cast recording 

***CTRL+D*** or ***CTRL+C*** to stop recording
***ENTER*** or ***CTRL+C*** to stop recording

![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/asciinema_stop.png)
<!-- blank line -->
----
<!-- blank line -->

### 8. exit docker container

```
exit
```

![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/asciinema_exit.png)
<!-- blank line -->
----
<!-- blank line -->

### 9. example screen-cast

[![asciicast](https://asciinema.org/a/453642.svg)](https://asciinema.org/a/453642?t=5)
