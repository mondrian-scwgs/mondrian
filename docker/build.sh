#!/bin/bash

TYPE=$1
PUSH=$2
CACHE=$3
USERNAME=$4
PASSWORD=$5


VERSION=`git describe --tags $(git rev-list --tags --max-count=1)`

cd $TYPE


CMD="docker build --build-arg VERSION=$VERSION -t quay.io/mondrianscwgs/$TYPE:$VERSION ."


if [ $CACHE == "Y" ]
then
    $CMD
else
    $CMD --no-cache
fi


if [ $PUSH == "Y" ]
then
    docker login quay.io -u $USERNAME --password $PASSWORD
    docker push quay.io/mondrianscwgs/$TYPE:$VERSION
fi


cd ../