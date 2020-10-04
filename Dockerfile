# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## create directories
RUN mkdir -p /01_data
RUN mkdir -p /02_code
RUN mkdir -p /03_output

RUN apt-get update \ 
	&& DEBIAN_FRONTEND=noninteractive apt-get install -y \
		libxml2-dev \
		libssl-dev \
		libcurl4-openssl-dev \
		libssh2-1-dev \
		libglfw3-dev \
		libglu1-mesa-dev \
		mesa-common-dev \
		freeglut3-dev \
		libfreetype6-dev \
		xorg \
		libx11-dev
	
## copy files
COPY /02_code/install_packages.R /02_code/install_packages.R
COPY /02_code/Fishplot_Docker.R /02_code/Fishplot_Docker.R

## install R-packages
RUN Rscript /02_code/install_packages.R

## run fishplot script
RUN Rscript /02_code/Fishplot_Docker.R