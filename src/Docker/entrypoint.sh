#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback

if [ -n "$LOCAL_USER_ID" ] # User / docker
then
    echo "User / docker"
    export DISPLAY=$DISPLAY
    USER_ID=${LOCAL_USER_ID:-9001}
    echo "Starting with User : $USER"
    useradd $USER --shell /bin/bash -u $USER_ID -o -c "" -g sudo -m  
    export HOME=/home/$USER
    echo "$USER ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
    exec /usr/sbin/gosu $USER /app/Diagnostic/DeCo.py $DATA_DIR $PARAM_FILE
elif [[ `whoami` == "root" ]] # Root / docker
then
    Xvfb :1 -screen 0 1024x768x16 &
    export DISPLAY=$DISPLAY
    echo "Root / docker"
    echo "Starting with UID : `whoami`"  # Root / docker
    exec "/app/../DeCo.py $DATA_DIR $PARAM_FILE"
else
    echo "User / singularity"
    echo "Starting with UID : `whoami`"  # User / singularity
    exec "/app/DeCo.py $DATA_DIR $PARAM_FILE"
fi


