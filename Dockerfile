# Our clusters run Ubuntu and so shall we
FROM ubuntu:16.04
LABEL maintainer="University of Michigan Center for Statistical Genetics"
WORKDIR /code
COPY . /code


#### Dependency installation (several options provided as examples- uncomment the one you need)
### Default GCC version that comes with the container
RUN apt-get update && \
    apt-get -y install apt-utils make gcc g++ gfortran zlib1g-dev r-mathlib

# Installing a newer GCC version requires adding a PPA.
# Provide a symlink for cases where current makefiles hardcode `g++`
#RUN apt-get update && apt-get install -y  apt-utils software-properties-common && \
#    add-apt-repository ppa:ubuntu-toolchain-r/test && \
#    apt-get update && \
#    apt-get -y install make gcc-6 g++-6 gfortran zlib1g-dev r-mathlib && \
#    ln -s `readlink /usr/bin/g++-6` /usr/bin/g++

## Analysis workflow tools. Only needed if you want to actually run the program in the container, eg to verify results
# RUN apt-get -y install vim tabix

## Build the executable. If this fails, the container will fail to build.
RUN cd libStatGen && make clean && cd ../libRareMetal && make clean && cd .. && \
    make

# For now just start a shell and keep it running, so user can connect to the shell inside the container for their work.
CMD /bin/bash