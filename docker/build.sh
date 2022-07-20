#!/bin/bash
set -e

TYPE=$1
USERNAME=$2
PASSWORD=$3
BETA=$4

VERSION=`git describe --tags $(git rev-list --tags --max-count=1)`

if [ $BETA == "Y" ]
then 
    DOCKER_VERSION="${VERSION}beta"
else
    DOCKER_VERSION=$VERSION
fi

cd $TYPE


docker build --build-arg VERSION=$VERSION -t quay.io/mondrianscwgs/$TYPE:$DOCKER_VERSION .



docker login quay.io -u $USERNAME --password $PASSWORD
docker push quay.io/mondrianscwgs/$TYPE:$DOCKER_VERSION


cd ../
