FROM continuumio/miniconda3
RUN mkdir /code
WORKDIR /code

# ADD requirements.txt /code/
# RUN conda install --file requirements.txt
RUN conda install django
RUN conda install -c bioconda bwa
RUN conda install -c conda-forge biopython
RUN conda install -c anaconda pandas
RUN conda install -c anaconda numpy
RUN conda install -c anaconda gunicorn

ADD . /code/
