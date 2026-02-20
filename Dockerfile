FROM debian:sid-slim

WORKDIR /app/

RUN chsh -s /bin/bash

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    rpm \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2025.12-2-Linux-aarch64.sh -O Anaconda.sh 
RUN /bin/bash Anaconda.sh -b -p /opt/conda 
RUN rm -f Anaconda.sh 
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh 
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc 
RUN echo "conda activate base" >> ~/.bashrc 
RUN /opt/conda/bin/conda clean -afy 

ENV PATH="/opt/conda/bin:/opt/conda/sbin:${PATH}"

RUN pip install gitdir==1.2.7

#RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh"

RUN /opt/conda/bin/conda init bash
RUN /opt/conda/bin/conda tos accept

RUN wget --quiet https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-aarch64-linux.tar.gz -O ncbi-blast-2.17.0+-aarch64-linux.tar.gz
RUN tar -xzf ncbi-blast-2.17.0+-aarch64-linux.tar.gz -C ./
RUN cp ncbi-blast-2.17.0+/bin/* /usr/local/bin/

RUN /opt/conda/bin/conda install -c bioconda mafft -y 
RUN /opt/conda/bin/conda install -c bioconda nextclade=3.18.1 -y
RUN /opt/conda/bin/conda install -c conda-forge biopython=1.86 -y
RUN /opt/conda/bin/conda install -c conda-forge pandas=3.0.0 -y
RUN /opt/conda/bin/conda install -c conda-forge networkx=3.6.1 -y

RUN pip install virallc

CMD [ "virallc" ]