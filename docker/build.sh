TYPE=$1
VERSION=$2


cd $TYPE

docker build --build-arg VERSION=$VERSION -t quay.io/mondrianscwgs/$TYPE:$VERSION . --no-cache
docker push quay.io/mondrianscwgs/$TYPE:$VERSION

cd ../