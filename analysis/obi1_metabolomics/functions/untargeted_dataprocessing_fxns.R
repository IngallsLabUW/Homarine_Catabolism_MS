#
# Copyright: Copyright 2022, Integral Consulting Inc. All rights reserved.
#
# Purpose:
#
# Project Information:
#   Name:
#   Number:
#
# History:
# Date	     Remarks
# YYYY-MM-DD
#

# PACKAGES & SPECIAL FUNCTIONS ----
library(tidyverse)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    ID_1 = rownames(cormat)[row(cormat)[ut]],
    ID_2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


# Define function to read in MS-DIAL Data----
#mode should be the analytical fraction ("HILICPos", "HILICNeg", or "RP")

MSDIAL_read <- function(file1, Mode) {
  dat_out <- read_delim(file1,
             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4,
             show_col_types = FALSE,
             name_repair = "unique_quiet")%>%
    mutate(Column = Mode) %>%
    mutate("ID" = as.character(`Alignment ID`)) %>%
    rename("RT" = `Average Rt(min)`,
           'mz' = `Average Mz`,
           "MS2_bool" = `MS/MS assigned`,
           "MS2" = `MS/MS spectrum`,
           "Name" = `Metabolite name`,
           "Adduct" = `Adduct type`,
           "Note" = `Post curation result`,
           "Fill" = `Fill %`,
           "Ref_RT" = `Reference RT`,
           "Ref_mz" = `Reference m/z`,
           "RT_matched" = `RT matched`,
           "mz_matched" = `m/z matched`,
           "MS2_matched" = `MS/MS matched`,
           "SN_ave" = `S/N average`) %>%
    select(ID, everything()) %>%
    select(-`Alignment ID`)
  return(dat_out)
}


# Adduct, isotope, and pos/neg match finder function-----
### adduct finding function
find_adducts <- function(ID_num, dataset, adduct_list, RT_ad_tol, adduct.error.ppm) {
  dataset.2 <- dataset %>%
    mutate(merge.key = "x")

  MF.limits <- dataset %>%
    filter(ID == ID_num) %>%
    mutate(RT_low = RT - RT_ad_tol,
           RT_high = RT + RT_ad_tol)
  ad.limits <- cbind(MF.limits, adduct_list)%>%
    mutate(adduct.mass = ((mz) - (1.007276*polarity))/abs(Charge) + Mass.change,
           adduct.mass.high = adduct.mass+(adduct.error.ppm/10^6)*adduct.mass,
           adduct.mass.low = adduct.mass-(adduct.error.ppm/10^6)*adduct.mass) %>%
    mutate(merge.key = "x") %>%
    select(merge.key, Ion, RT_low, RT_high, adduct.mass.low, adduct.mass.high, adduct.mass, ad.polarity) %>%
    rename("polarity" = ad.polarity)
  ad.detect <- full_join(dataset.2, ad.limits, by = c("polarity", "merge.key"))%>%
    select(!merge.key)%>%
    rowwise() %>%
    filter(RT >= RT_low & RT <= RT_high)%>%
    filter(mz >= adduct.mass.low & mz <= adduct.mass.high) %>%
    rename("adduct.ID" = ID) %>%
    mutate(ID = ID_num) %>%
    mutate(ad.ppm = abs((adduct.mass-mz)/adduct.mass*10^6))
  return(ad.detect)
}


# WRAPPER Dereplication function ------

dereplicate_MFs <- function(
    dat_filename1,
    dat_type_1,
    dat_filename2 = NA,
    dat_type_2 = NA,
    comment_tag = "m",
    sample_tag = "2112",
    drop_tags = c("DDA", "_Blk_"),
    RT_tol = 0.2, #minutes to allow for possible adduct/isotope
    cor_tol = 0.95,  # minimum correlation coefficient for considering if adduct/iso
    add_ppm_tol = 20, # adduct mass tolerance (in ppm)
    adduct_file = "raw_input/KHeal_adduct_list.csv",
    date_tag = ""
) {
  # GET POLARITIES
  if (str_detect(dat_type_1, "Pos")) {
    polarity_1 = 1
  } else {
    polarity_1 = -1
  }
  if (str_detect(dat_type_2, "Pos")) {
    polarity_2 = 1
  } else {
    polarity_2 = -1
  }


  # LOAD DATA   ----
  dat_adducts <- read_csv(adduct_file, show_col_types = FALSE) %>%
    rename(ad.polarity = polarity)
  dat_1 <- MSDIAL_read(dat_filename1, dat_type_1) %>%
    # keep only good peaks
    filter(Comment == comment_tag) %>%
    mutate(polarity = polarity_1)
  if (!is.na(dat_filename2)) {
    dat_2 <- MSDIAL_read(dat_filename2, dat_type_2) %>%
      # keep only good peaks
      filter(Comment == comment_tag) %>%
      mutate(polarity = polarity_2)
  }

  # FIND CORRELATED MFs ----
  ## Combine dat1 and dat2 ----
  if (is.na(dat_filename2)) {
    dat_cmb <- dat_1%>%
      mutate(ID = paste0(Column, "_", ID))
  } else {
    dat_cmb <- dat_1 %>%
      bind_rows(dat_2) %>%
      mutate(ID = paste0(Column, "_", ID))
  }

  ## Make matrix of areas ----
  # remove blanks and DDA
  dat_m <- dat_cmb %>%
    select(ID, contains(sample_tag)) %>%
    select(-c(matches(paste(drop_tags, collapse = "|"))))%>%
    as.matrix()
  row.names(dat_m) <- dat_m[, 1]
  dat_m <- dat_m[, -1]
  ## Perform and clean up correlation results ----
  corr_results <- rcorr(t(dat_m), type = c("pearson"))
  dat_corr <- flattenCorrMatrix(corr_results$r, corr_results$P) %>%
    left_join(dat_cmb %>%
                select(ID, RT, mz, polarity) %>%
                rename(
                  RT_1 = RT,
                  mz_1 = mz,
                  polarity_1 = polarity
                ),
              by = c("ID_1" = "ID")
    ) %>%
    left_join(dat_cmb %>%
                select(ID, RT, mz, polarity) %>%
                rename(
                  RT_2 = RT,
                  mz_2 = mz,
                  polarity_2 = polarity
                ),
              by = c("ID_2" = "ID")
    ) %>%
    mutate(
      RT_diff = abs(RT_1 - RT_2),
      mz_diff = abs(mz_1 - mz_2)
    )
  ### Filter out only the very correlated mass features ----
  dat_vcorr <- dat_corr %>%
    filter(
      abs(RT_diff) < RT_tol,
      cor > cor_tol
    )

  ## For each individual MFs, see if they are adducts or isotopes-----
  list_vcorr <- dat_vcorr %>% split(dat_vcorr$ID_1)

  # Loop around list (or map, eventually)
  found_adducts <- tibble()
  message(paste0("Searching over ", length(list_vcorr), " groups of MFs"))
  for (i in 1:length(list_vcorr)) {
    comp.ID <- unique(list_vcorr[[i]]$ID_1)
    ID.dat <- list_vcorr[[i]] %>%
      select(ID_2, mz_2, RT_2, polarity_2)  %>%
      rename(
        ID = ID_2,
        mz = mz_2,
        RT = RT_2,
        polarity = polarity_2
      ) %>%
      distinct() %>%
      bind_rows(list_vcorr[[i]] %>%
                  select(ID_1, mz_1, RT_1, polarity_1) %>%
                  rename(
                    ID = ID_1,
                    mz = mz_1,
                    RT = RT_1,
                    polarity = polarity_1
                  ) %>%
                  filter(polarity == unique(list_vcorr[[i]]$polarity_1)) %>%
                  distinct())
    # this function asks for each compound - if this is M+H or M-H, do we see any matching adducts or isotopes?
    adduct.out <- find_adducts(
      comp.ID,
      ID.dat,
      dat_adducts,
      RT_ad_tol = RT_tol,
      adduct.error.ppm = add_ppm_tol
    )
    if (length(adduct.out$adduct.ID) > 1) {
      found_adducts <- bind_rows(found_adducts, adduct.out)
#      print(paste0(as.character(i), " found an adduct!"))
    }
  }
  # CLEAN UP ANNOATIONS ------
  ## Get summary of IDs that are adducts ------
  found_adducts2 <- found_adducts %>%
    mutate(ad.ppm = abs(ad.ppm)) %>%
    group_by(adduct.ID, ID) %>%
    filter(ad.ppm == min(ad.ppm)) %>%
    select(ID, adduct.ID, mz, RT, ad.ppm, Ion) %>%
    distinct()

  adduct_summary <- found_adducts2 %>%
    filter(ad.ppm != 0) %>%
    group_by(adduct.ID) %>%
    filter(ad.ppm == min(ad.ppm)) %>%
    rename(
      psued_ID = ID,
      ID = adduct.ID
    ) %>%
    mutate(add_annotation = paste0(Ion, " of ", psued_ID))

  adduct_summary2 <- adduct_summary %>%
    group_by(ID) %>%
    summarise(add_annotation = paste0(add_annotation, collapse = " or "))

  ## Get summary of IDs that are psuedo molecular ions ------
  pseudos <- found_adducts2 %>%
    filter(ad.ppm == 0) %>%
    filter(ID %in% adduct_summary$psued_ID) %>%
    mutate(add_annotation = Ion) %>%
    ungroup() %>%
    select(ID, add_annotation)


  ## Add to dataframe -----
  dat_cmb2 <- dat_cmb %>%
    left_join(adduct_summary2, by = "ID") %>%
    rows_update(pseudos, by = "ID") %>%
    select(ID, RT, mz, MS2, add_annotation, starts_with(date_tag))


  return(dat_cmb2)
}




# BMIS (MSDial input) ------
BMIS_MSdialoutput <-
  function(dat1_filename,
           sample_key,
           is_names_filename,
           is_dat1_filename,
           is_dat2_filename = NA,
           cut_off1,
           cut_off2,
           smps_to_dump = c(),
           is_to_dump = c()) {
    # Import data -----
    dat1 <- read_csv(dat1_filename,
                     show_col_types = FALSE)
    sample_key_all <- read_csv(sample_key, show_col_types = FALSE)


    # pivot MSdial data ----
    dat <- dat1 %>%
      pivot_longer(
        -c(ID, RT, mz, MS2, add_annotation),
        names_to = "replicate_name", values_to = "area"
      ) %>%
      filter(!is.na(area)) %>%
      rename(mass_feature = ID)
    # drop blanks, standards, bad IS, and bad samps
    dat <- dat %>%
      filter(!str_detect(replicate_name, "_Blk_")) %>%
      filter(!str_detect(replicate_name, "_Std_")) %>%
      mutate(area = as.numeric(area)) %>%
      filter(!str_detect(add_annotation, " of ") | is.na(add_annotation)) %>%
      filter(!(replicate_name %in% smps_to_dump))

    # prep is data ----
    is_names <- read_csv(
      is_names_filename,
      show_col_types = FALSE)
    if (is.na(is_dat2_filename)){
      is_dat <- read_csv(
        is_dat1_filename,
        show_col_types = FALSE
      )
    } else {
      is_dat <- read_csv(
        is_dat1_filename,
        show_col_types = FALSE
      ) %>%
        bind_rows(read_csv(is_dat2_filename,
                           show_col_types = FALSE
        ))
    }
    browser()
    is_dat <- is_dat %>%
      clean_names()    %>%
      rename(mass_feature = precursor_ion_name) %>%
      filter(mass_feature %in% is_names$is_name)    %>%
      filter(replicate_name %in% dat$replicate_name)    %>%
      mutate(area = as.numeric(area)) %>%
      suppressWarnings() %>%
      select(replicate_name, mass_feature, area) %>%
      filter(!(mass_feature %in% is_to_dump))

    samp_key_is <- sample_key_all %>%
      filter(replicate_name %in% is_dat$replicate_name) %>%
      select(replicate_name, is_fraction) %>%
      filter(!is.na(is_fraction)) %>%
      mutate(
        mass_feature = "Inj_vol",
        area = is_fraction,
        replicate_name = replicate_name
      ) %>%
      select(replicate_name, area, mass_feature)

    is_dat <- bind_rows(is_dat, samp_key_is)

    # Look at extraction replication of the Internal Standards----
    IS_inspectPlot <-
      ggplot(is_dat, aes(x = replicate_name, y = area)) +
      geom_bar(stat = "identity") +
      facet_wrap(~mass_feature, scales = "free_y") +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 5
        ),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)
      ) +
      ggtitle("IS Raw Areas")

    dat <- dat %>%
      mutate(QC_area = area) %>% #TODO - add in some QC
      select(replicate_name, mass_feature, area, QC_area) %>%
      mutate(date = str_extract(replicate_name, "^\\d*"))
    is_dat <- is_dat %>%
      mutate(date = str_extract(replicate_name, "^\\d*"))

    # Calculate mean values for each IS----
    is_means <- is_dat %>%
      group_by(mass_feature, date) %>%
      summarise(ave = mean(as.numeric(area)), .groups = "drop") %>%
      mutate(ave = ifelse(mass_feature == "Inj_vol", 1, ave))


    # Normalize to each internal Standard----
    split_dat <- list()
    browser()
    for (i in 1:length(unique(is_dat$mass_feature))) {
      split_dat[[i]] <-
        bind_rows(is_dat %>%
                    mutate(QC_area = area), dat) %>%
        mutate(MIS = unique(is_dat$mass_feature)[i]) %>%
        left_join(
          is_dat %>%
            rename(MIS = mass_feature, is_area = area) %>%
            select(MIS, replicate_name, is_area),
          by = c("replicate_name", "MIS")
        ) %>%
        left_join(is_means %>%
                    rename(MIS = mass_feature),
                  by = c("date", "MIS")
        ) %>%
        mutate(adjusted_area = QC_area / is_area * ave)
    }
    dat_norm <- bind_rows(split_dat) %>% select(-is_area, -ave)

    # Break Up the Names (Name structure must be:  Date_type_ID_replicate_anythingextraOK)----
    dat_norm <- dat_norm %>%
      separate(
        replicate_name,
        into = c(
          "runDate",
          "type", "samp_id", "replicate"
        ),
        sep = "_",
        remove = FALSE,
        extra = "drop",
        fill = "right"
      )


    # Find the B-MIS for each mass_feature----
    # Look only the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
    # then choose which IS reduces the RSD the most (Poo.Picked.IS)
    poodat <- dat_norm %>%
      filter(type == "Poo") %>%
      group_by(samp_id, mass_feature, MIS) %>%
      summarise(RSD_ofPoo_IND = sd(adjusted_area,
                                   na.rm = TRUE
      ) / mean(adjusted_area, na.rm = TRUE), .groups = "drop") %>%
      mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
      group_by(mass_feature, MIS) %>%
      summarise(RSD_ofPoo = mean(RSD_ofPoo_IND, na.rm = TRUE))

    poodat <- poodat %>% left_join(poodat %>%
                                     group_by(mass_feature) %>%
                                     summarise(
                                       poo_picked_is =
                                         unique(MIS)[which.min(RSD_ofPoo)][1]
                                     ),
                                   by = "mass_feature"
    )

    # Get the starting point of the RSD (Orig_RSD), calculate the change in the RSD, say if the MIS is acceptable----
    poodat <- poodat %>%
      left_join(
        poodat %>%
          filter(MIS == "Inj_vol") %>%
          mutate(Orig_RSD = RSD_ofPoo) %>%
          select(-RSD_ofPoo, -MIS),
        by = c("mass_feature", "poo_picked_is")
      ) %>%
      mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
      mutate(percent_change = del_RSD / Orig_RSD) %>%
      mutate(accept_MIS = (percent_change > cut_off1 &
                             Orig_RSD > cut_off2))



# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----
# Adds a column that has the BMIS, not just Poo.picked.IS
# Changes the finalBMIS to inject_volume if its no good
fixedpoodat <- poodat %>%
  filter(MIS == poo_picked_is) %>%
  mutate(FinalBMIS = ifelse((accept_MIS == "FALSE"), "Inj_vol", poo_picked_is))
newpoodat <-
  poodat %>%
  left_join(fixedpoodat %>% select(mass_feature, FinalBMIS), by = c("mass_feature")) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)
Try <- newpoodat %>% filter(FinalBMIS != "Inj_vol")
QuickReport <- paste0(
  round(length(Try$mass_feature) / length(newpoodat$mass_feature), digits = 3)*100,
  " % of MFs  picked a BMIS. RSD improvement cutoff = ",
  cut_off1,
  ".  RSD minimum cutoff = ",
  cut_off2
)

# Evaluate the results of your BMIS cutoff-----
IS_toISdat <- dat_norm %>%
  filter(mass_feature %in% is_dat$mass_feature) %>%
  select(mass_feature, MIS, adjusted_area, type) %>%
  filter(type == "Smp") %>%
  group_by(mass_feature, MIS) %>%
  summarise(RSD_ofSmp = sd(adjusted_area) / mean(adjusted_area)) %>%
  left_join(
    poodat %>% select(mass_feature, MIS, RSD_ofPoo, accept_MIS),
    by = c("mass_feature", "MIS")
  )

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol")


ISTest_plot <- ggplot() +
  geom_point(
    dat = IS_toISdat,
    shape = 21,
    color = "black",
    size = 2,
    aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS)
  ) +
  scale_fill_manual(values = c("white", "dark gray")) +
  geom_point(
    dat = injectONlY_toPlot,
    aes(x = RSD_ofPoo, y = RSD_ofSmp),
    size = 3
  ) +
  facet_wrap(~mass_feature) +
  theme_bw()

# Get all the data back - and keep only the MF-MIS match set for the BMIS----
# Add a column to the longdat that has important information from the FullDat_fixed,
# then only return data that is normalized via B-MIS normalization
BMIS_normalizedData <-
  newpoodat %>%
  select(mass_feature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(dat_norm %>% rename(FinalBMIS = MIS), by = c("mass_feature", "FinalBMIS")) %>%
  unique() %>%
  filter(!mass_feature %in% is_dat$mass_feature)

BMISlist <-
  list(
    IS_inspectPlot,
    QuickReport,
    ISTest_plot,
    BMIS_normalizedData
  )

return(BMISlist)
  }
