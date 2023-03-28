# latch base image + dependencies for latch SDK --- removing these will break the workflow
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9c8f-main
run pip install latch==2.14.2
run mkdir /opt/latch

# install R requirements
run apt-get update --yes && \
    apt-get install --yes software-properties-common && \
    add-apt-repository "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev wget
copy environment.R /opt/latch/environment.R
run Rscript /opt/latch/environment.R

# copy all code from package (use .dockerignore to skip files)
copy . /root/

# latch internal tagging system + expected root directory --- changing these lines will break the workflow
arg tag
env FLYTE_INTERNAL_IMAGE $tag
workdir /root
