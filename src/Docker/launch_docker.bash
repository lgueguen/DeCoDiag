#!/bin/bash

##
# Script used to launch docker run
#

##
# Launch with 1st argument as the working directory and 2nd argument
# as parameter file as seen from the working directory

# if 2nd argument is void, window interface is opened to ask

# if 1st argument is void, actual directory is used as working
# directory

if [ -z "$1" ] # if no argument
then
    DATA_DIR=`eval echo ~$USER`
    ATTACH_DIR=home
else
    DATA_DIR=`readlink -f $1`
    ATTACH_DIR=data
fi

docker run -it -e DISPLAY=$DISPLAY -e USER=$USER -e LOCAL_USER_ID=`id -u $USER` \
       -v /tmp/.X11-unix:/tmp/.X11-unix -v $DATA_DIR:/app/$ATTACH_DIR \
       -e PARAM_FILE=$2 -e DATA_DIR=/app/$ATTACH_DIR \
       --cpus 1 \
       --rm --name DeCo anal_dock


