FROM --platform=linux/amd64 nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update \
    && DEBIAN_FRONTEND=noninteractive \
    && apt-get install --no-install-recommends -y \
    git \
    wget \
    curl

# RUN apt install ca-certificates

WORKDIR /home/repos

# # Install ColabFold
# RUN wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
# RUN nvcc --version
# RUN gcc --version
# RUN bash install_colabbatch_linux.sh
# RUN echo "export PATH=\"/home/repos/localcolabfold/colabfold-conda/bin:$PATH\"" >> ~/.bashrc

# # Install FragFold

RUN apt install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt install python3.7 python3-pip git -y
RUN pip3 install --upgrade pip

# # Final version, once the repo is public
# RUN apt install ssh -y
# RUN git clone -v git@github.com:swanss/FragFold.git \
#   && cd FragFold

# Temporary workaround
WORKDIR FragFold 
COPY . .
SHELL ["/bin/bash", "-c"]
RUN python3 -m pip install -e .

# ENV PATH /root/miniconda3/bin:$PATH

# RUN mkdir -p ~/miniconda3 \
#     && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh \
#     && bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 \
#     && rm -rf ~/miniconda3/miniconda.sh \
#     && ~/miniconda3/bin/conda init bash \
#     && . /root/miniconda3/etc/profile.d/conda.sh \
#     && conda create --name fragfold -y \
#     && conda activate fragfold

# RUN pip3 install -e .