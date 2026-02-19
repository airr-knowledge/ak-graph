FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install package dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-dev \
        python3-pip \
        python3-sphinx \
        python3.12-venv \
        build-essential \
        libomp-dev \
        cmake \
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

# compile C++ code
RUN mkdir /ak-graph
COPY . /ak-graph
RUN cd /ak-graph/compairr && make install
RUN cd /ak-graph && make

# Create and enable virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install networkit inside venv
RUN pip install --upgrade pip
RUN pip install networkit matplotlib numpy pandas powerlaw networkx seaborn scipy scikit-learn psycopg

RUN mkdir /work
WORKDIR /work
