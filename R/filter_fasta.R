# filter_fasta

# The minimal_taxa setting in Y&S Step 6 "Paralogy pruning" scripts filter paralog trees
# by a minimal number of taxa without considering ingroup / outgroup status (except for
# RT). Use this script to further filter the results from filter_1to1_orthologs.py,
# prune_paralogs_MI.py etc by outgroup/ingroup status.

# given a set of fasta files with some minimum TOTAL taxa cutoff from an ortholog pruning step,
# filter the fasta files down to only those including a minimum number of INGROUP taxa only

library(ape)
library(xlsx)

setwd("/Volumes/Storage/thelypteris/working/yang_smith/run7")

# function to extract cluster name(s) from alignment, including "cluster" part and "1rr" part
extract_cluster_name <- function (align) {
  # grab either names or rownames of alignment depending if its in matrix format or a list
  if (length(dim(align)) == 0) {
    otu <- names(align)
  } else if (length(dim(align)) > 0) {
    otu <- rownames(align)
  }

  otu <- otu[grep("cluster", otu)]
  otu <- sapply(otu, function (x) strsplit(x, split="_"))
  otu_first <- sapply(otu, function (x) x[2])
  otu_second <- sapply(otu, function (x) x[3])
  otu <- paste(otu_first, otu_second, sep="_")
  cluster_name <- unique(otu)
  return(cluster_name)
}

# get associated data for 1kp samples
onekp_data <- read.xlsx("/Volumes/Transcend/Documents/Projects/JSPS/one_kp/1kP-Sample-List.xlsx",
                        sheetIndex = 1,
                        colIndex = 1:19,
                        startRow=1,
                        endRow=84,
                        stringsAsFactors=FALSE)

ingroup <- onekp_data$Code[onekp_data$out_in == "IN"]
ingroup <- ingroup[!(is.na(ingroup))]

onekp_data$Family_subfam <- paste(onekp_data$Family, onekp_data$subfamily, sep="_" )
onekp_data$Family_subfam <- gsub("_NA", "", onekp_data$Family_subfam)

onekp_data$genus <- sapply(onekp_data$Species, function (x) strsplit(x, split="_")[[1]][1])

eu2_family_list <- sort(unique(onekp_data$Family[onekp_data$out_in == "IN"]))
eu2_subfamily_list <- sort(unique(onekp_data$Family_subfam[onekp_data$out_in == "IN"]))
eu2_genus_list <- sort(unique(onekp_data$genus[onekp_data$out_in == "IN"]))


filter_fasta <- function (MCL_settings, homolog_type, taxaset, filter_type="min_taxa", min_taxa, xinping_seqs) {
  fasta.names <- list.files(paste(MCL_settings, "/", homolog_type, "/fasta/", sep=""))
  fasta.names <- fasta.names[grep(".fa$", fasta.names)]
  fasta.files <- lapply(paste(MCL_settings, "/", homolog_type, "/fasta/", fasta.names, sep=""), read.dna, format="fasta")

  # filter fasta files based on having at least one species per eupoly II family with seq length within 1 sd of overall mean seq length
  # input: candidate bait fasta file
  # output: T/F whether the candidate passes the filter
  filter_eupoly2_family <- function (alignment) {
    # minimum cutoff seq length is mean overall seq length - 1 sd
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just EuII taxa
    eu2_alignment <- alignment[names(alignment) %in% ingroup]
    # subsetted alignment with just EuII taxa having min seq length
    eu2_alignment <- eu2_alignment[sapply(eu2_alignment, length) > min_length ]
    # get list of families for subsetted alignment
    eu2_families <- onekp_data$Family[match(names(eu2_alignment), onekp_data$Code)]
    # check if all 9 eu2 families are included at least once in subsetted alignment
    all(eu2_family_list %in% eu2_families)
  }

  # filter fasta files based on having at least one species per eupoly II subfamily with seq length within 1 sd of overall mean seq length
  # input: candidate bait fasta file
  # output: T/F whether the candidate passes the filter
  filter_eupoly2_subfamily <- function (alignment) {
    # minimum cutoff seq length is mean overall seq length - 1 sd
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just EuII taxa
    # eu2_alignment <- alignment[full_subfamilies %in% eu2_subfamily_list]
    eu2_alignment <- alignment[names(alignment) %in% ingroup]
    # subsetted alignment with just EuII taxa having min seq length
    eu2_alignment <- eu2_alignment[sapply(eu2_alignment, length) > min_length ]
    # get list of families for subsetted alignment
    eu2_subfamilies <- onekp_data$Family_subfam[match(names(eu2_alignment), onekp_data$Code)]
    # check if all 12 eu2 subfamilies are included at least once in subsetted alignment
    all(eu2_subfamily_list %in% eu2_subfamilies)
  }

  # filter fasta files based on having at least one species per eupoly II subfamily with seq length within 1 sd of overall mean seq length
  # input: candidate bait fasta file
  # output: T/F whether the candidate passes the filter
  filter_eupoly2_genus <- function (alignment) {
    # get list of genera for each seq in alignment
    # full_genera <- onekp_data$genus[match(names(alignment), onekp_data$Code)]
    # minimum cutoff seq length is mean overall seq length - 1 sd
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just EuII taxa
    #eu2_alignment <- alignment[full_genera %in% eu2_genus_list]
    eu2_alignment <- alignment[names(alignment) %in% ingroup]
    # subsetted alignment with just EuII taxa having min seq length
    eu2_alignment <- eu2_alignment[sapply(eu2_alignment, length) > min_length ]
    # get list of genera for subsetted alignment
    eu2_genera <- onekp_data$genus[match(names(eu2_alignment), onekp_data$Code)]
    # check if all 24 eu2 subfamilies are included at least once in subsetted alignment
    all(eu2_genus_list %in% eu2_genera)
    # length(which(eu2_genus_list %in% eu2_genera))
  }

  filter_xinping <- function (alignment_names, xinping_seqs) {
    clusters <- sapply(alignment_names, function (x) strsplit(x, split="_"))
    clusters_first <- sapply(clusters, function (x) x[1])
    clusters_second <- sapply(clusters, function (x) x[2])
    clusters <- paste(clusters_first, clusters_second, sep="_")
    clusters %in% xinping_seqs
  }

  if (filter_type == "min_taxa" && !is.na(min_taxa)) {
    # filter fasta files based on minimum number of taxa in ingroup
    count.in <- lapply(fasta.files, FUN = function (x) length(names(x)[names(x) %in% taxaset]))
    count.total <- lapply(fasta.files, FUN = function (x) length(names(x)))
    pass_filter <- count.in >= min_taxa
  } else if (filter_type == "min_taxa" && is.na(min_taxa)) {
    stop("need to provide minimum number of taxa for filtering")
  } else if (filter_type == "family") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_family)
  } else if (filter_type == "subfamily") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_subfamily)
  } else if (filter_type == "genus") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_genus)
  } else if (filter_type == "xinping") {
    pass_filter <- filter_xinping (fasta.names, xinping_seqs)
  } else {
    stop ("need to choose valid filtering method")
  }

  fasta.filtered <- fasta.files[pass_filter]
  fasta.names.filtered <- fasta.names[pass_filter]

  results <- list (
    fasta.filtered = fasta.filtered,
    fasta.names.filtered = fasta.names.filtered,
    fasta.filtered.length = length(fasta.filtered),
    fasta.raw.length = length(fasta.files)
  )
  return(results)
}

# first run check_yang_smith (check taxon occupancy: across ingroup taxa only) to look at taxon occupancy chart and choose best min taxa cutoff

# run7: for comparison, run for each filtering scheme (min taxa, family, and subfamily)
# the "filtered_fasta" folder will contain all of the fasta files, but I can subselect using the filtered lists
MCL_settings <- "hit-frac0.3_I1.4_e5"
prune_settings <- "ortho_121"
homolog_type <- "ortho_121"

# 121 pruning, hit-frac0.3_I1.4_e5, only those that are also in xinping 159 genes
# read names of clusters that are also in xinping 159 gene set (obtained by running compare_baits.R)
setwd("/Volumes/Storage/thelypteris/working/yang_smith/compare_baits/common_baits")
clusters_in_xinping <- list.files()
clusters_in_xinping <- lapply(clusters_in_xinping, read.dna, format="fasta")
clusters_in_xinping <- lapply(clusters_in_xinping, extract_cluster_name)
# some alignments have multiple clusters. unlist these.
clusters_in_xinping <- unlist(clusters_in_xinping)
setwd("/Volumes/Storage/thelypteris/working/yang_smith/run7")

filtered_fasta_results <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", taxaset=ingroup, filter="xinping", xinping_seqs=clusters_in_xinping)
filtered_fasta_results$fasta.filtered.length
filtered_fasta_results$fasta.raw.length
# write list of filtered cluster names for loop scripts
write(filtered_fasta_results$fasta.names.filtered, file=(paste0(MCL_settings, "/", prune_settings, "/introns/filtered_align_list_xinping")))

# 121 pruning, hit-frac0.3_I1.4_e5, min taxa = 26, should result in 1660 orthogroups
# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run7")
filtered_fasta_results <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", taxaset=ingroup, filter="min_taxa", min_taxa=26)
filtered_fasta_results$fasta.filtered.length
filtered_fasta_results$fasta.raw.length
# write list of filtered cluster names for loop scripts
write(filtered_fasta_results$fasta.names.filtered, file=(paste0(MCL_settings, "/", prune_settings, "/introns/filtered_align_list_mintaxa")))
# write out filtered fasta files
# setwd(paste0("/Volumes/Storage/thelypteris/working/yang_smith/run7/", MCL_settings, "/", prune_settings, "/filtered_fasta"))
# mapply(write.fastas, filtered_fasta_results$fasta.filtered, filtered_fasta_results$fasta.names.filtered)

# at least one taxon per subfamily
filtered_fasta_results <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", filter="subfamily")
filtered_fasta_results$fasta.filtered.length
# write list of filtered cluster names for loop scripts
write(filtered_fasta_results$fasta.names.filtered, file=(paste0(MCL_settings, "/", prune_settings, "/introns/filtered_align_list_subfamily")))

# at least one taxon per family
filtered_fasta_results <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", filter="family")
filtered_fasta_results$fasta.filtered.length
# write list of filtered cluster names for loop scripts
write(filtered_fasta_results$fasta.names.filtered, file=(paste0(MCL_settings, "/", prune_settings, "/introns/filtered_align_list_family")))

# at least one taxon per genus - ONLY 3 ALIGNMENTS!
filtered_fasta_results <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", filter="genus")
filtered_fasta_results$fasta.filtered.length
filtered_fasta_results$fasta.raw.length
# write list of filtered cluster names for loop scripts
# write(filtered_fasta_results$fasta.names.filtered, file=(paste0(MCL_settings, "/", prune_settings, "/introns/filtered_align_list_mintaxa")))






# (old versions)

# try a more conservative cutoff (want about 500 genes total)
# 121 pruning, hit-frac0.3_I1.4_e5, min taxa = 28, should result in 903 orthogroups
#one2one_hitfrac0.3_I1.4_cut28 <- filter_fasta("hit-frac0.3_I1.4_e5", "ortho_121", taxaset=ingroup, min_taxa=28)
#one2one_hitfrac0.3_I1.4_cut28$fasta.filtered.length
#one2one_hitfrac0.3_I1.4_cut28$fasta.raw.length

# 121 pruning, hit-frac0.4_I2_e5, min taxa = 27, should result in 1829 orthogroups
# (intron filtering will probably cut this number in half)
# one2one_hitfrac0.4_I2_cut27 <- filter_fasta("hit-frac0.4_I2_e5", "ortho_121", taxaset=ingroup, min_taxa=27)
# one2one_hitfrac0.4_I2_cut27$fasta.filtered.length
# one2one_hitfrac0.4_I2_cut27$fasta.raw.length
#
# # 121 pruning, hit-frac0.3_I2_e5, min taxa = 26, should result in  1841 orthogroups
# # (intron filtering will probably cut this number in half)
# MCL_settings <- "hit-frac0.3_I2_e5"
# prune_settings <- "ortho_121"
# filtered_fasta_results <- filter_fasta(MCL_settings, prune_settings, taxaset=ingroup, min_taxa=26)
# filtered_fasta_results$fasta.filtered.length
# filtered_fasta_results$fasta.raw.length
#
# # 121 pruning, hit-frac0.4_I1.4_e5, min taxa = 27 -> 1542 orthogroups
# # (intron filtering will probably cut this number in half)
# MCL_settings <- "hit-frac0.4_I1.4_e5"
# prune_settings <- "ortho_121"
# filtered_fasta_results <- filter_fasta(MCL_settings, prune_settings, taxaset=ingroup, min_taxa=27)
# filtered_fasta_results$fasta.filtered.length
# filtered_fasta_results$fasta.raw.length




# # write list of filtered cluster names for loop scripts
# write(MO_cut10$fasta.names.filtered, file="6_introns/MO_align_list")
#
# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run5/ortho_MO/filtered_fasta")
# mapply(write.fastas, MO_towrite$fasta.filtered, MO_towrite$fasta.names.filtered)
#
# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run5/ortho_RT/filtered_fasta")
# mapply(write.fastas, RT$fasta.filtered, RT$fasta.names.filtered)

# run6
# write list of filtered cluster names for loop scripts
# write(MO_cut25$fasta.names.filtered, file="5_introns/MO_align_list")
# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run6/hit-frac0.4_I1.4_e5/ortho_MO/filtered_fasta")
# mapply(write.fastas, MO_cut25$fasta.filtered, MO_cut25$fasta.names.filtered)



# run6
# MO method, min taxa=25 results in ca. 2800 orthogroups
# MO_cut25 <- filter_fasta("hit-frac0.4_I1.4_e5", "ortho_MO", taxaset=ingroup, min_taxa=25)
# MO_cut25$fasta.filtered.length
# MO_cut25$fasta.raw.length

# run5
# # MO method, min taxa=10 results in ca. 2000 orthogroups
# MO_cut10 <- filter_fasta("ortho_MO", taxaset=ingroup, min_taxa=10)
# MO_cut10$fasta.filtered.length
# MO_cut10$fasta.raw.length
#
# # MO method, min taxa=11 results in ca. 1000 orthogroups
# MO_cut11 <- filter_fasta("ortho_MO", taxaset=ingroup, min_taxa=11)
# MO_cut11$fasta.filtered.length
# MO_cut11$fasta.raw.length

# I want to add the min taxa = 10 fasta files
# but I already wrote out min_taxa = 11 fasta files
# to save on time for aligning, only write out those in min taxa = 10 that aren't already in 11
# MO_towrite <- list()
# MO_towrite$fasta.filtered <- MO_cut10$fasta.filtered[!(MO_cut10$fasta.names.filtered %in% MO_cut11$fasta.names.filtered)]
# MO_towrite$fasta.names.filtered <- MO_cut10$fasta.names.filtered[!(MO_cut10$fasta.names.filtered %in% MO_cut11$fasta.names.filtered)]
#
# # RT, min taxa = 11
# RT <- filter_fasta("ortho_RT", taxaset=ingroup, min_taxa=11)
# RT$fasta.filtered.length
# RT$fasta.raw.length



#write(one2one_hitfrac0.3_I1.4_cut26$fasta.names.filtered, file="hit-frac0.3_I1.4_e5/ortho_121/introns/one2one_hitfrac0.3_I1.4_cut26_align_list")
#write(one2one_hitfrac0.4_I2_cut27$fasta.names.filtered, file="hit-frac0.4_I2_e5/ortho_121/introns/one2one_hitfrac0.4_I2_cut27_align_list")

# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run7/hit-frac0.3_I1.4_e5/ortho_121/filtered_fasta")
# mapply(write.fastas, one2one_hitfrac0.3_I1.4_cut26$fasta.filtered, one2one_hitfrac0.3_I1.4_cut26$fasta.names.filtered)

# setwd("/Volumes/Storage/thelypteris/working/yang_smith/run7/hit-frac0.4_I2_e5/ortho_121/filtered_fasta")
# mapply(write.fastas, one2one_hitfrac0.4_I2_cut27$fasta.filtered, one2one_hitfrac0.4_I2_cut27$fasta.names.filtered)
