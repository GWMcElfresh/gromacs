FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

#dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    wget \
    git \
    libnetcdf-dev \
    libopenmpi-dev openmpi-bin \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

#setup conda
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

#set conda environment variables
ENV PATH="/opt/conda/bin:$PATH"
RUN conda init && conda config --set always_yes yes

#activate conda environment (python 3.9)
RUN conda create -n openmm-env python=3.9 && \
    echo "conda activate openmm-env" >> ~/.bashrc

#get latest ambertools from conda
RUN conda install -n openmm-env -c conda-forge openmm ambertools mdtraj numpy scipy pandas matplotlib seaborn

#verify install:
RUN conda run -n openmm-env python -c "import openmm; print(openmm.version.version)"

#reset to default workspace
WORKDIR /workspace
CMD ["/bin/bash"]
