FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install package dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-sphinx \
        python3-scipy \
        python3-networkx \
	    python3-pandas \
	    python3-matplotlib \
        python3-igraph \
        python3-cairocffi \
        build-essential \
        software-properties-common \
        autoconf \
        automake \
        libtool \
        pkg-config \
        ca-certificates \
        libssl-dev \
        wget \
        git \
        curl \
        language-pack-en \
        locales \
        locales-all \
        nano \
        gdb \
        valgrind \
        libboost-all-dev \
        libpqxx-dev libpq-dev
