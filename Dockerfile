FROM rocker/shiny:4.3.3

# System dependencies for R packages (magick, curl, xml2, tidyverse graphics stack)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libmagick++-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages — single layer for cache efficiency.
# install.packages() only warns on failure, so verify afterwards and fail the
# build if any package is missing (e.g. tidyverse silently dropping ragg).
RUN R -e "pkgs <- c( \
    'shiny', 'visNetwork', 'plotly', 'igraph', 'reactable', 'DT', \
    'readxl', 'writexl', 'tidyverse', 'scales', 'shinyWidgets', \
    'shinyjs', 'htmlwidgets', 'magick'); \
  install.packages(pkgs, repos='https://cran.r-project.org/', Ncpus=parallel::detectCores()); \
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
  if (length(missing)) stop('Failed to install: ', paste(missing, collapse=', ')); \
  cat('All R packages installed successfully\n')"

# heatmaply drives the clustered gene-list heatmap in the Transcriptome Browser.
# Kept in its own layer: appending it to the list above would invalidate that
# layer and recompile tidyverse from source on every rebuild.
RUN R -e "install.packages('heatmaply', repos='https://cran.r-project.org/', Ncpus=parallel::detectCores()); \
  if (!'heatmaply' %in% rownames(installed.packages())) stop('Failed to install: heatmaply'); \
  cat('heatmaply installed successfully\n')"

# App is mounted as a volume at runtime — no COPY needed
EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838, launch.browser=FALSE)"]
