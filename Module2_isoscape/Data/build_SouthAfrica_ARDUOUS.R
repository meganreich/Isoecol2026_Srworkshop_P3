###############################################################################
# Build the harmonised and screened bioavailable 87Sr/86Sr dataset for South
# Africa (ARDUOUS format) used in the ISOECOL 2026 isoscape workshop.
#
# Inputs : Data/SrBaseline.csv (global, ARDUOUS format)
#          Data/combined_other_bioavailable_sr_data_clean.csv (Wang compilation
#          + IsoBank)
#          Project_Raster/clipped_za/r.srsrq1.tif (study-area mask)
# Output : Data/SouthAfrica_bioavailable_Sr_ARDUOUS.csv
#
# Screening rules
#  1. Records with coordinates and 87Sr/86Sr, within the South Africa rasters,
#     0.700 < 87Sr/86Sr < 0.800
#  2. Bioavailable proxies only: plants, soils, and modern fauna with small home
#     ranges (rodents, hyraxes, hares, elephant shrews, snails, tortoises,
#     insects, small territorial antelopes, small carnivores)
#  3. Removed: waters (catchment-integrated), rocks/sediments, fish, humans and
#     hominins, primates, elephants/ivory, domestic cattle, birds, large and
#     migratory herbivores, large carnivores, fossil specimens (Lehmann et al.
#     2018 Elandsfontein; Balter et al. 2012 Sterkfontein/Swartkrans; fossil
#     rodents of Copeland et al. 2010), unidentified fauna
#  4. Duplicates (same measurement in both files) removed, keeping the
#     SrBaseline record; when SrBaseline coordinates are rounded to 0.01 degree,
#     the more precise coordinates of the Wang compilation are used
###############################################################################

library(terra)

base  <- read.csv("Data/SrBaseline.csv", check.names = FALSE)
other <- read.csv("Data/combined_other_bioavailable_sr_data_clean.csv", check.names = FALSE)
arduous_cols <- names(base)

# ---- 1. Bring the "other" compilation into the ARDUOUS columns --------------
names(other) <- sub("^X87Sr86Sr", "87Sr86Sr", names(other))
other$scientific_name <- trimws(gsub(" ", " ", other$scientific_name))

# Free-text material descriptions -> ARDUOUS material_type vocabulary
mat_map <- c("Plant"        = "organism : plant : plant tissue",
             "Tooth enamel" = "organism : animal : animal tissue : tooth : enamel",
             "Teeth"        = "organism : animal : animal tissue : tooth",
             "Dentin"       = "organism : animal : animal tissue : tooth : dentin",
             "Bone"         = "organism : animal : animal tissue : bone",
             "bone"         = "organism : animal : animal tissue : bone",
             "Bone/Ivory"   = "organism : animal : animal tissue : tusk",
             "Shell"        = "organism : animal : animal tissue : shell",
             "Snail shell"  = "organism : animal : animal tissue : shell",
             "Tortoise"     = "organism : animal : animal tissue",
             "Quill"        = "organism : animal : animal tissue",
             "Insect"       = "organism : animal : whole animal")
free <- other$material_type %in% names(mat_map)
other$material_type[free] <- mat_map[other$material_type[free]]
other$analysis_type[is.na(other$analysis_type) | other$analysis_type == ""] <- "bulk solid"
# give the Wang records an identifier
no_id <- is.na(other$sample_measurement_id) | other$sample_measurement_id == ""
other$sample_measurement_id[no_id] <- sprintf("WANG-%04d", seq_len(sum(no_id)))

for (col in setdiff(arduous_cols, names(other))) other[[col]] <- NA
other <- other[, arduous_cols]

# Combine: SrBaseline first so that its (richer) records are kept for duplicates
sr <- rbind(base, other)
sr$`87Sr86Sr` <- suppressWarnings(as.numeric(sr$`87Sr86Sr`))

# ---- 2. Coordinates, value range, study area ----------------------------------
sr <- sr[!is.na(sr$collection_decimal_latitude) & !is.na(sr$collection_decimal_longitude) &
         !is.na(sr$`87Sr86Sr`) & sr$`87Sr86Sr` > 0.700 & sr$`87Sr86Sr` < 0.800, ]
mask <- rast("Project_Raster/clipped_za/r.srsrq1.tif")
pts  <- project(vect(sr, geom = c("collection_decimal_longitude", "collection_decimal_latitude"),
                     crs = "EPSG:4326"), mask)
sr <- sr[!is.na(extract(mask, pts, ID = FALSE)[, 1]), ]

# ---- 3. Bioavailable proxies -------------------------------------------------
mat   <- sr$material_type
taxon <- tolower(sr$scientific_name)
cit   <- sr$related_publication_citation

is_plant <- grepl("^organism : plant", mat)
is_soil  <- grepl("soil", mat)
is_fauna <- grepl("^organism : animal", mat)

small_range_taxa <- c("rodentia", "rhabdomys", "rhabodomys", "otomys", "gerbilliscus", "cryptomys",
                      "micaelamys", "pedetes", "xerus", "hystrix", "procavia", "lepus", "leptus",
                      "pronolag", "leporidae", "macroscel", "achatina", "metachatina", "cochlitoma",
                      "tropidophora", "agate snail", "bark snail", "genusnata", "gastropoda",
                      "chersina", "psammobates", "stigmochelys", "testudines", "tortoise",
                      "raphicerus", "sylvicapra", "pelea", "tragelaphus scriptus", "suricata",
                      "otocyon")
small_fauna <- is_fauna & Reduce(`|`, lapply(small_range_taxa, grepl, x = taxon, fixed = TRUE)) &
               !grepl("fossil", taxon)

# fossil assemblages (any material) and wide-ranging animals
fossil_study <- grepl("Lehmann", cit) | grepl("Balter", cit)

keep <- (is_plant | is_soil | small_fauna) & !fossil_study

# duplicates of excluded records (e.g. fossil rodents flagged only in the Wang
# compilation) are excluded too
key <- paste(round(sr$collection_decimal_latitude, 3), round(sr$collection_decimal_longitude, 3),
             round(sr$`87Sr86Sr`, 4))
keep <- keep & !(key %in% key[!keep & grepl("fossil", taxon)])

sr  <- sr[keep, ]
key <- key[keep]

# ---- 4. Duplicates -----------------------------------------------------------
# a) exact duplicates (same coordinates and 87Sr/86Sr)
dup <- duplicated(key)
sr  <- sr[!dup, ]

# b) the same measurement in both files, but SrBaseline coordinates rounded to
#    0.01 degree: same 87Sr/86Sr (to 1e-5) and coordinates within 0.011 degree.
#    We keep the SrBaseline record and give it the more precise coordinates.
is_wang <- grepl("^WANG", sr$sample_measurement_id)
drop <- rep(FALSE, nrow(sr))
used <- rep(FALSE, nrow(sr))
for (i in which(is_wang)) {
  j <- which(!is_wang & !used &
             abs(sr$`87Sr86Sr` - sr$`87Sr86Sr`[i]) < 1e-5 &
             abs(sr$collection_decimal_latitude - sr$collection_decimal_latitude[i]) <= 0.011 &
             abs(sr$collection_decimal_longitude - sr$collection_decimal_longitude[i]) <= 0.011)
  if (length(j) > 0) {
    j <- j[1]
    sr$collection_decimal_latitude[j]  <- sr$collection_decimal_latitude[i]
    sr$collection_decimal_longitude[j] <- sr$collection_decimal_longitude[i]
    used[j] <- TRUE
    drop[i] <- TRUE
  }
}
cat("Near-duplicates removed:", sum(drop), "\n")
sr <- sr[!drop, ]

write.csv(sr, "Data/SouthAfrica_bioavailable_Sr_ARDUOUS.csv", row.names = FALSE, na = "")

cat(nrow(sr), "records\n")
print(table(ifelse(grepl("^organism : plant", sr$material_type), "plant",
              ifelse(grepl("soil", sr$material_type), "soil", "fauna"))))
