FROM debian:bookworm-20241202


# Install R and default to using Posit package manager
# See: https://packagemanager.posit.co/client/#/repos/cran/setup
RUN apt-get update \
    && apt-get install -y gpg \
    && gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
    && gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' > /etc/apt/trusted.gpg.d/cran_debian_key.asc \
    && echo "deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/" > /etc/apt/sources.list.d/r.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev \
    && echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/bookworm/2024-12-16"))' >> /etc/R/Rprofile.site \
    && echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' >> /etc/R/Rprofile.site


# Basic apt installs
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        time \
        cmake \
        curl \
        zstd \
        lz4 \
        tabix \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgsl-dev \
        libharfbuzz-dev \
        libhdf5-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
        libglpk-dev

# Install R packages
RUN Rscript -e 'install.packages("pak")' \
    && echo \
            arrow \
            bench \
            BiocManager \
            devtools \
            hdf5r \
            RcppHNSW \
            remotes \
            RSpectra \
            Seurat \
            tidyverse \
            tiledb \
        | MAKEFLAGS="-j 4" Rscript -e 'pak::pkg_install(strsplit(readLines("stdin"), " ")[[1L]])' \
    && echo \
            chromVAR \
            DelayedArray \
            DelayedMatrixStats \
            HDF5Array \
            motifmatchr \
            BSgenome.Hsapiens.UCSC.hg19 \
            BSgenome.Hsapiens.UCSC.hg38 \
        | MAKEFLAGS="-j 4" Rscript -e 'BiocManager::install(strsplit(readLines("stdin"), " ")[[1L]], version="3.20")' \
    && echo \
            immunogenomics/presto@7636b3d0465c468c35853f82f1717d3a64b3c8f6 \
            mojaveazure/seurat-disk@877d4e18ab38c686f5db54f8cd290274ccdbe295 \   
        | MAKEFLAGS="-j 4" Rscript -e 'remotes::install_github(strsplit(readLines("stdin"), " ")[[1L]], repos=)' \
    && rm -rf /root/.cache/R /tmp/*

# Install ArchR library
RUN Rscript -e 'remotes::install_github("GreenleafLab/ArchR@d9e741c980c7c64e5348c97a74d146cc95f8ba76", repos = BiocManager::repositories(version="3.20"))' \
    && rm -rf /root/.cache/R /tmp/*

# Install tiledbsoma library
RUN apt-get install -y --no-install-recommends git \
    && MAKEFLAGS="-j 4" Rscript -e "install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))" \
    && MAKEFLAGS="-j 4" Rscript -e "install.packages('cellxgene.census', repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))" \
    && rm -rf /root/.cache/R /tmp/*

# Install micromamba + scanpy  + other python dependencies
COPY conda-env.yml requirements.txt /root

ENV ENV_NAME="bpcells"
ENV MAMBA_ROOT_PREFIX="/root/micromamba"
ENV MAMBA_EXE="/root/.local/bin/micromamba"
RUN curl -L micro.mamba.pm/install.sh | bash \
    && bash -c '. "/root/micromamba/etc/profile.d/mamba.sh" \
              && micromamba env create -f /root/conda-env.yml \
              && micromamba activate $ENV_NAME \
              && pip install -r /root/requirements.txt'

# Install lzbench
RUN curl -L https://github.com/inikep/lzbench/archive/7534ba679d07c2594e4a710c4dab3fbcea80ee41.zip > lzbench.zip \
    && unzip lzbench.zip \
    && rm lzbench.zip \
    && mv lzbench-7534ba679d07c2594e4a710c4dab3fbcea80ee41 lzbench-master \
    && cd lzbench-master \
    && make -j 4 \
    && cd / 

# Install BPCells at the v0.3.0 release
RUN MAKEFLAGS="-j 4" Rscript -e 'remotes::install_github("bnprks/BPCells",  subdir="r", ref="e7551359f85b20ca448481a1ba913ea947b6966e"); quit(status=!requireNamespace("BPCells"))' \
    && rm -rf /root/.cache/R /tmp/*

# Make it so the conda env is activated by default
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
