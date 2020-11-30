FROM rocker/shiny:3.5.2

# System libraries of general use and some required for the R packages
RUN apt-get update && apt-get install -y \
  libssl-dev \
  libxml2-dev

# System library dependency for the genetic interaction tool app
RUN R -e "install.packages(c('devtools'))"
RUN R -e "install.packages('dplyr', repo='https://mran.microsoft.com/snapshot/2020-06-24/')"
RUN R -e "install.packages('shiny', repo='https://mran.microsoft.com/snapshot/2020-06-24/')"
RUN R -e "install.packages('cgdsr', repo='https://mran.microsoft.com/snapshot/2020-10-24/')"
RUN R -e "install.packages('mclust', repo='https://mran.microsoft.com/snapshot/2020-02-24/')"
RUN R -e "install.packages('stringr', repo='https://mran.microsoft.com/snapshot/2020-10-24/')"
RUN R -e "install.packages('XML', repo='https://mran.microsoft.com/snapshot/2020-06-24/')"
RUN R -e "install.packages('seqinr', repo='https://mran.microsoft.com/snapshot/2020-08-24/')"
RUN R -e "install.packages('tidyr', repo='https://mran.microsoft.com/snapshot/2020-06-24/')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "require(BiocManager)"
RUN R -e "BiocManager::install(version='3.8')"
RUN R -e "BiocManager::install('data.table')"
RUN R -e "BiocManager::install('maftools')"

# Copy the app to the image
WORKDIR /srv/shiny-server
COPY setup.R .
COPY data data/
COPY www www/
COPY functions.R .
COPY app.R .

# Setup from R
RUN Rscript setup.R

# Allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server

# Run app
CMD ["/usr/bin/shiny-server.sh"]
