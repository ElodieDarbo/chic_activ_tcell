# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.0.3

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# copy necessary files
## app folder
#COPY /scripts ./scripts
COPY /renv_0.13.2.tar.gz ./renv_0.13.2.tar.gz

## data folder
#COPY /data ./data

RUN mkdir /data

VOLUME ["/data"]

RUN mkdir /scripts

VOLUME ["/scripts"]


## renv.lock file
COPY /renv.lock ./renv.lock

# install renv & restore packages
RUN Rscript -e 'install.packages("/renv_0.13.2.tar.gz",repos=NULL,type="source")'
RUN Rscript -e 'renv::restore()'

RUN rm /renv_0.13.2.tar.gz

# expose port
EXPOSE 3838
# run app on container start
CMD ["echo","The app is available at http:/localhost:3838"]
CMD ["R", "-e", "shiny::runApp('/scripts', host = '0.0.0.0', port = 3838)"]