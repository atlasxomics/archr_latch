FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

# Install R
RUN apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" && \
    apt-get install -y \
        r-base \
        r-base-dev \
        apt-transport-https \
        build-essential \
        gfortran \
        libatlas-base-dev \
        libbz2-dev \        
        libcurl4-openssl-dev \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgit2-dev \
        libgsl-dev \
        libicu-dev \
        liblzma-dev \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libpcre3-dev \
        libssl-dev \
        libtcl8.6 \
        libtiff5 \
        libtk8.6 \
        libxml2-dev \
        libxt-dev \
        libx11-dev \
        locales \
        make \
        pandoc \
        tzdata \
        zlib1g-dev        

# Have to install devtools like this; see https://stackoverflow.com/questions/20923209, also cairo
RUN apt-get install -y r-cran-devtools

# Installation of R packages with renv
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/renv_1.0.5.tar.gz', repos = NULL, type = 'source')"
COPY renv.lock /root/renv.lock
COPY .Rprofile /root/.Rprofile
RUN mkdir /root/renv
COPY renv/activate.R /root/renv/activate.R
COPY renv/settings.json /root/renv/settings.json
RUN R -e "renv::restore()"

# Install python packages
COPY requirements.txt /root/requirements.txt
RUN python3 -m pip install -r requirements.txt

COPY custom_ArchR_genomes_and_annotations/ /root/custom_ArchR_genomes_and_annotations/

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
