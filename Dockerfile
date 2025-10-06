FROM continuumio/miniconda3

# Install blaze
RUN pip3 install blaze2 

# Install minimap2 and samtools
RUN apt-get update && apt-get install -y minimap2 samtools && rm -rf /var/lib/apt/lists/*

# install R + bambu
# Install R via conda 
RUN conda install -y -c conda-forge r-base r-essentials && conda clean -afy

# Install bambu from Bioconductor
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
          BiocManager::install('bambu', ask=FALSE, update=FALSE)"

RUN curl --proto '=https' --tlsv1.2 -LsSf https://github.com/COMBINE-lab/oarfish/releases/latest/download/oarfish-installer.sh | sh 
RUN cp /root/.cargo/bin/oarfish /usr/local/bin/ && chmod +x /usr/local/bin/oarfish

CMD ["bash"]

