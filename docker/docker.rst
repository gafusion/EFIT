Docker Tools
===============

Using pre-built Docker containers registered at gitlab
------------------------------------------------------

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

Running a container locally
---------------------------

It is wholly possible to run pipeline job locally on your computer
using gitlab's runner instead of in the cloud or using one of the
efit shared runners. Gitlab provides a runner specificaly for this
purpose, called "gitlab-runner". Installation instructions for
gitlab-runner can be found here (different platforms are indicated
on the left side of the web page):

https://docs.gitlab.com/runner/install/osx.html

Alternatively, one can install the docker image and run the container,
on the local system volume mount, which is simpler and fast. See

https://docs.gitlab.com/runner/install/docker.html

for more details. From the above page, 

.. code:: bash

    docker run -d --name gitlab-runner --restart always \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /srv/gitlab-runner/config:/etc/gitlab-runner \
    gitlab/gitlab-runner:latest

On MacOS, replace "srv" above with "/Users/Shared". Once the container is
running in docker, you can run a pipeline job locally with the following
command:

.. code:: bash

    gitlab-runner exec docker <job_name>

where <job_name> is the name of the individual job in the pipeline.
"gitlab-runner" should be installed in "/usr/local/bin/gitlab-runner".


Building and registering docker containers with gitlab
------------------------------------------------------

Basic steps for registering at gitlab:

.. code:: bash

    docker login registry.gitlab.com   # Only need to do once.  Use gitlab creds
    DNAME=aocc
    docker build -f Dockerfile.$DNAME -t registry.gitlab.com/efit-ai/efit/$DNAME .
    docker push registry.gitlab.com/efit-ai/efit/$DNAME

Alternatively to build and run it, you can use `docker-compose`::

    docker-compose up --build -d efitai-$DNAME  # build the container using Dockerfile and start it

Alternative to docker compose
-----------------------------


To start a container from the EFIT pre-built Docker image on GitLab
(eg. clang below), run:

.. code:: bash

    cd docker
    docker pull registry.gitlab.com/efit-ai/efit/clang
    docker images   # Show all images in your local docker

With that image available, you can start a container with an interactive shell:

.. code:: bash

    docker run -ti registry.gitlab.com/efit-ai/efit/clang /bin/bash


Note that this doesn't have the directories mounted which is why we use
docker-compose.


Creating your own docker image
------------------------------

To create your own image and push it to the efit-ai project, edit a file
called "Dockerfile.<image_name>" where "<image_name>" is the name of the
image that you want to create, eg. "clang". You can use one of the existing
docker image files located in "efit-ai/efit/docker" as a template. These
files will soon be moved to "efit-ai/docker".

.. code:: bash

     docker build -f Dockerfile.<image_name> -t registry.gitlab.com/efit-ai/efit/<image_name> .
     docker images
     docker login registry.gitlab.com 
     docker push registry.gitlab.com/efit-ai/efit/<image_name>

It is important when building the image to tag it ("-t") as is done above. This
will obviate the need to tag the image when pushing to the registry. Once built,
the image should be mounted in docker. One must login to the registry prior to
pushing the image.

Note that if at the end of the push, you see the error:
"unauthorized: authentication required" then the push was unsuccessful. This seems
to be caused by a timeout after logging in to the registry. Since images can be
very large, this can happen if you have a slower upload connection (10 MBps is
often not enough).

You can check the container registry at the URL at the top of this document. If you
see your new image there with no tags, then the push was not successful. If there is
one tag, then you are good to go.

