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
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Install packages
RUN R -e "install.packages(c('Cairo', 'BiocManager', 'Matrix', 'Seurat'))"
RUN R -e "devtools::install_github('immunogenomics/harmony')"
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R -e "library('ArchR'); ArchR::installExtraPackages()"

RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"

RUN python3 -m pip install slims-python-api

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN echo sup
RUN python3 -m pip install --upgrade latch

# DO NOT REMOVE
run pip install latch==2.36.3
run mkdir /opt/latch

# Copy workflow data (use .dockerignore to skip files)
copy . .latch/* /root/

# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root