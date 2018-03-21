# Using the Docker image
A Docker image is provided for development purposes, in order to ensure that the package can be built from a 
clean environment. This *is not* intended for production use of Raremetal.

It may also be useful for, eg, running on unsupported operating systems such as Mac OS. 
(subject to any resource limitations of the host VM)

If developing on Mac OS, you may need to install [Docker for Mac](https://docs.docker.com/docker-for-mac/). 

## Commands

### Build the image
`docker build . -t raremetal2img`

### Run the image
The image may be converted into a running container instance as follows:
`docker run -it -d --name raremetal2container raremetal2img`

Because the default command is to run a bash shell, it will remain running and can be addressed by name for further 
commands. Eg, you can access an interactive shell to try the executable via:
`docker exec -it raremetal2container  /bin/bash`

Or `docker stop raremetal2container` the container, then `docker start raremetal2container` later.


### Clean up
To clean up:

Remove unused/ stopped containers:
`docker container prune`

Remove unused build images (eg when rebuilding):
`docker image prune`

