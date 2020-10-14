# fishplot-docker
### This repository contains a dockerfile to build a docker container for an interegraged pipeline of [FishPlot](https://github.com/chrisamiller/fishplot), [ClonEvol](https://github.com/hdng/clonevol) and [SciClone](https://github.com/genome/sciclone) packages. 

The base image used for this container is rocker/r-base. Here are the steps to build and run the docker container.

## Folder Structure

The folders are as follow -
- 01_data : This folder contains the data file
- 02_code : This folder contains the code
- 03_output : The output of the pipeline will be generated in this folder

## Get the fishplot-docker

1) Clone the repository `git clone https://github.com/genespectrum/fishplot-docker.git`
2) Navigate to the folder fishplot-docker

## Building and running the image

3) Build the image using the base image from Docker Hub. The below command will build the container on a local installation.  Navigate to the fishplot-docker folder and then run this command. It will take some time to build the container. 

       `docker build -t fishdocker .`

4) Once finished, run the container using the below command. 

       `docker run -it --rm -v /path/of/01_data:/01_data -v /path/of/03_output:/03_output fishdocker`

Using the `-v` argument signals Docker which local folders to map to the created folders inside the container. This is important because we want to both get our dataframe inside the container and save our output from the workflow locally so it isn’t lost once the container is stopped.

The output files generated from this pipeline can be accessed from the 03_output folder. 
