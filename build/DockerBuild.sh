#!/bin/bash


docker build -t wenjin27/deeplc:gpu ./
docker push wenjin27/deeplc:gpu
docker run --user 1000:1000 wenjin27/deeplc:gpu



