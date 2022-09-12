Docker Tools
===============

Using pre-built Docker containers registered at gitlab
-------------------------------------------------------------

We have pre-built docker containers.  They are used for both the CI pipeline
("shared" runners which run in the cloud) and for ease in development.  For
development, they enable developers to checkout an image, launch a bash shell,
and start developing in a common environment shared with all developers.  The
instructions below assume that you have docker installed. 

Here we show how to checkout the docker from the registry, launch it, and then
start a bash session into it.  We will use `clang` as build system, but note
that there are other images in the registry as seen here:

https://gitlab.com/efit-ai/efit/container_registry


Also, we use `docker-compose` to simplify some of the commands and make a nicer
environment.  We include the raw docker commands at the bottom for contrast.
To start a container from the EFIT pre-built Docker image on GitLab, run:

.. code:: bash

    cd docker
    echo COMPOSE_PROJECT_NAME=$USER > .env  # [optional] specify a project name
    DNAME=aocc
    docker-compose pull efitai-$DNAME   # Get the most updated image from the repo
    docker images   # Show all images in your local docker.  Should have efit
    docker-compose up -d efitai-$DNAME  # start the container 
    docker ps                           # Show the running container

Then to launch an interactive Bash session inside that container, do:

.. code:: bash

    docker-compose exec efitai-$DNAME bash

By default, you are placed into `/efit` directory, which is mounted from your
local file system.  This allows you to the edit on your local file system, but
still develop in your local container.  To build a version of efit in your
container, do:
  
.. code:: bash

    mkdir ../build   # Need to be in container file system
    cd ../build
    ../efit/docker/config_docker_$DNAME.sh
    make
    make test

This should build and pass tests since this is tested nightly.

After exiting the container, do not forget to cleanup after yourself:

.. code:: bash

    exit
    docker-compose down # stop and remove the container

Building and registering docker containers with gitlab
-------------------------------------------------------------

Creating your own docker image::

    docker-compose  TODO

Basic steps for registering at gitlab::

    docker login registry.gitlab.com   # Only need to do once.  Use gitlab creds
    DNAME=aocc
    docker build -f Dockerfile.$DNAME -t registry.gitlab.com/efit-ai/efit/$DNAME .
    docker push registry.gitlab.com/efit-ai/efit/$DNAME

Alternatively to build and run it, you can use `docker-compose`::

    docker-compose up --build -d efitai-$DNAME  # build the container using Dockerfile and start it

Alternative to docker compose
-----------------------------


To start a container from the EFIT pre-built Docker image on GitLab, run:

.. code:: bash

    cd docker
    docker pull registry.gitlab.com/efit-ai/efit/clang
    docker images   # Show all images in your local docker

With that image available, you can start a container with an interactive shell:

.. code:: bash

    docker run -ti registry.gitlab.com/efit-ai/efit/clang /bin/bash


Note that this doesn't have the directories mounted which is why we use
docker-compose
