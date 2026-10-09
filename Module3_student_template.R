# =============================================================================
#  IsoEcol 2026 Strontium isotope geolocation for terrestrial ecology: Module 3
#  isotope-based geographic assignment
# =============================================================================
#
#  Online tutorial: https://meganreich.github.io/Isoecol2026_Srworkshop_P3/
#
#  How to use this script:
#   1. Follow the tutorial in your browser.
#   2. Copy each code block (use the copy button in its top-right corner)
#      and paste it under the matching section below.
#   3. Run lines with Ctrl + Enter (Cmd + Enter on a Mac).
#   4. Add your own notes with # so you can come back to this later.


# 0. Setup ---------------------------------------------------------------------
# Load the packages we'll use. Install any that are missing with
# install.packages("package_name")
install.packages(c("tidyverse", "assignR","terra","isocat","rnaturalearth", "sf", "tidyterra", "viridis"))

library(tidyverse)
library(assignR)
library(terra)
library(rnaturalearth)
library(sf)
library(tidyterra)
library(viridis)
library(isocat)

