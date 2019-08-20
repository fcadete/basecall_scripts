FROM ubuntu:18.04

#RUN apt-get update \
#    && apt-get install -y wget lsb-release \
#    && export PLATFORM=$(lsb_release -cs) \
#    && wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - \
#    && echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list \
#    && apt-get update \
#    && apt-get install -y ont-guppy-cpu

RUN apt-get update \
    && apt-get install -y wget lsb-release \
    && wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_2.3.7_linux64.tar.gz \
    && tar -xvf ont-guppy-cpu_2.3.7_linux64.tar.gz

RUN apt-get update \
    && apt-get install -y python3 python3-pip git build-essential ipython3 \
    && git clone https://github.com/rrwick/Porechop.git \
    && cd Porechop \
    && python3 setup.py install


RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base \
    && R -e "install.packages(c('data.table', 'futile.logger','ggplot2','optparse','plyr','readr','reshape2','scales','viridis','yaml'))" \
    && wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O MinIONQC.R


RUN apt-get update \
    && apt install -y default-jre \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
    && unzip fastqc_v0.11.8.zip \
    && chmod 755 FastQC/fastqc \
    && sed -i 's/kmer\s*ignore\s*1/kmer      ignore      0/g' /FastQC/Configuration/limits.txt

RUN git clone https://github.com/lh3/minimap2 \
    && cd minimap2 \
    && make

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure \
    && make \
    && make install

RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz \
    && tar -xzf BBMap_38.44.tar.gz

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz \
    && tar -xzf jellyfish-1.1.12.tar.gz \
    && cd jellyfish-1.1.12 \
    && ./configure \
    && make \
    && make install \
    && ldconfig

RUN git clone https://github.com/lh3/seqtk.git \
    && cd seqtk \
    && make

RUN apt-get update \
    && apt-get install -y curl libcurl4-gnutls-dev \
    && R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biostrings'); biocLite('Rsamtools')" \
    && R -e "install.packages('RcppRoll')"

RUN apt-get update \
    && apt-get install libxml2-dev \
    && R -e "install.packages('XML'); source('https://bioconductor.org/biocLite.R'); biocLite('rGADEM')"

RUN pip3 install beautifulsoup4

RUN apt-get install -y libssl-dev libgsl-dev \
    && R -e "source('https://bioconductor.org/biocLite.R'); biocLite('TFBSTools'); biocLite('JASPAR2014')"
    

RUN apt-get update \
    && apt-get install -y python-pip python2.7 \
    && pip install telomerehunter
 
