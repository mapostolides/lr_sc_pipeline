FROM continuumio/miniconda3

RUN pip3 install blaze2 

RUN apt-get update && apt-get install -y minimap2 samtools && rm -rf /var/lib/apt/lists/*

CMD ["bash"]

