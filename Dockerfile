FROM informaticsmatters/rdkit-python3-debian:Release_2019_09
ADD requirements.txt requirements.txt
USER root
RUN pip install -r requirements.txt
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y git procps
ADD setup.py README.rst /usr/local/fragmentor/
ADD frag /usr/local/fragmentor/frag
RUN pip install /usr/local/fragmentor
WORKDIR /usr/local/fragmentor/