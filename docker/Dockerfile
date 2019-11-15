FROM ubuntu:18.04

MAINTAINER dongqingsun96@gmail.com

WORKDIR /root/

ADD . /root/Annotation/

# install miniconda3
RUN apt-get update && apt-get install -y wget bzip2 git vim && \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
bash Miniconda3-latest-Linux-x86_64.sh -b && rm Miniconda3-latest-Linux-x86_64.sh

# redirct /bin/sh from dash to bash
RUN echo "dash dash/sh boolean false" | debconf-set-selections && DEBIAN_FRONTEND=noninteractive dpkg-reconfigure dash

ENV PATH=/root/miniconda3/bin:${PATH}

# install cellranger
# RUN mkdir /root/Software/ && cd /root/Software/ && wget -O cellranger-3.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz?Expires=1573619377&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU3MzYxOTM3N319fV19&Signature=enT9cujyI1NR3gDmFrGLvVWSZt-iA23rgwKYmbpFF4fem8CyG4TGBhuCW3D5QioWdpP071vGE~33Tgiw~XDd8n8~kGKizviLyu8YhH0YW9zWoT2MMv7gOnG82cLNhQ0BnaRgeXfMaz4cxGaPwmqCoD994Cq2iLjgaDC3tvqC3h7OfNpCuCpZxmhvIPj1i7UwiLpi8N~dY7pvHV2hDSV6xQtPpb2cOGOrxhnZPXvtiYMgAC8VbUaU8uetY7an4Inea~q2sPDLBQWd0jTwybROPSJZYrJfhCDWvfr~7xD~ZckyXzknMv-fNcZs5q-zknihjZ01fARDAfMQgwlsPjK7cw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
# wget -O cellranger-atac-1.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-1.1.0.tar.gz?Expires=1573619397&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWF0YWMvY2VsbHJhbmdlci1hdGFjLTEuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU3MzYxOTM5N319fV19&Signature=QtG56kghv~eAZ1R3Z6jT55dMgqAd~y9hCeMEdnPS7JgYw1KWy2L20Fx8Of5vpvbbuN8SGs4iBPoypqMLBxAyjD6zGhryozMbZ6U8YlscCY04unjAFo036~8hJi6~sGgMG~YsiwnQMO2RNjR-CE5ARF1I6GGuwksTM6SBzgPT9tgcdAZEcdSWJoQXsSltRLdKpbEcSasmsSXbFlMfpKT5l10Qi5mG94gj5BFgAdTpJDkdlvAl9C9RdbDRWgfNl9V6R0xqWVQt7lwEvMOc5jSWtwuqnIzkQQriigqRXZxhoB--KI7brG4sezGFVqdRov7op2lIpar0WygE1uo2VGvY1A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
# tar -xzvf cellranger-3.1.0.tar.gz && tar -xzvf cellranger-atac-1.1.0.tar.gz 

# download reference for cellranger
# RUN mkdir /root/Annotation/ && cd /root/Annotation/ && \
# wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz && \
# wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz && \
# wget http://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.1.0.tar.gz && \
# wget http://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.1.0.tar.gz && \
# tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz && \
# tar -xzvf refdata-cellranger-mm10-3.0.0.tar.gz && \ 
# tar -xzvf refdata-cellranger-atac-GRCh38-1.1.0.tar.gz && \
# tar -xzvf refdata-cellranger-atac-mm10-1.1.0.tar.gz 

# add conda channels and install softwares
RUN conda config --add channels defaults && conda config --add channels conda-forge && \
conda config --add channels bioconda && conda install -c liulab-dfci maestro

# install giggle and rabit
RUN apt install -y gcc make autoconf zlib1g-dev libbz2-dev libcurl4-openssl-dev libssl-dev ruby && \
mkdir /root/Software/ && cd /root/Software/ && git clone https://github.com/ryanlayer/giggle.git && cd giggle && make && \
echo "export PATH=/root/Software/giggle/bin:$PATH" >> /root/.bashrc && \
apt-get install -y build-essential libgsl-dev libblas-dev liblapack-dev && conda install -c conda-forge lapack && \
cd ../ && wget http://rabit.dfci.harvard.edu/download/Rabit.tar.gz && tar xvf Rabit.tar.gz && \
cd Rabit && ./configure --prefix=/root/miniconda3 LDFLAGS="-L/usr/local/lib -L/usr/lib -L/root/miniconda3/lib" CPPFLAGS="-I/usr/local/include -I/usr/include -I/root/miniconda3/include" && make && make check && make install

# download rabit index and giggle index
# RUN mkdir /root/Annotation/ && cd /root/Annotation/ && wget http://cistrome.org/~chenfei/MAESTRO/rabit.tar.gz && \
# wget http://cistrome.org/~chenfei/MAESTRO/giggle.tar.gz && \
RUN cd /root/Annotation && tar -xzvf rabit.tar.gz && tar -xzvf giggle.tar.gz && \
rm rabit.tar.gz && rm giggle.tar.gz

