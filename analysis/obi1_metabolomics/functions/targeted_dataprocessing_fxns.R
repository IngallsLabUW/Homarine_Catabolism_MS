#  targeted_dataprocessing_fxns.R
# _____________________________________________________________________________
# PURPOSE: -----
#  Example of a project database
#
#
# PROJECT INFORMATION:-----
#   Name:
#   Number:
#
# HISTORY:-----
# 	 Date		        Remarks
# _____________________________________________________________________________
# 	 20220420       Using this as a template: https://github.com/kheal/Gradients1_SemiTargeted3/blob/master/SourceCode/QC_HILICNeg_Cultures.R and adapting for these data
# _____________________________________________________________________________

library(tidyverse)
library(janitor)
library(here)
library(zoo)
library(snakecase)
library(RaMS)

# QC function ------------
QE_QC <- function(dat_filename,
                  std_flag = "Std",
                  blank_flag = "Blk",
                  fileout = paste0("QCd", dat_filename),
                  sn_min = 3,
                  ppm_flex = 6,
                  area_min = 40000,
                  rt_flex = 2.5,
                  blank_ratio_max = 0.3) {
  dat <- read_csv(dat_filename,
    show_col_types = FALSE
  )

  ## Clean up column types-----
  dat_og <- dat %>%
    clean_names() %>%
    select(-protein_name, -protein) %>%
    mutate_at(
      c(
        "retention_time",
        "area",
        "background",
        "height",
        "mass_error_ppm"
      ),
      as.numeric
    ) %>%
    suppressWarnings()

  ## Pull out samples, add sn, ppm, area min flags ---
  dat_flag <- dat_og %>%
    mutate(sn_flag = ifelse((area / background < sn_min), "sn_flag", NA)) %>%
    mutate(ppm_flag = ifelse((abs(mass_error_ppm) > ppm_flex), "ppm_flag", NA)) %>%
    mutate(areamin_flag = ifelse((area < area_min), "areamin_flag", NA))

  ## Get rts of standards----
  dat_rt_stds <- dat_og %>%
    filter(str_detect(replicate_name, std_flag)) %>%
    select(precursor_ion_name, retention_time) %>%
    group_by(precursor_ion_name) %>%
    summarise(rt_std = mean((retention_time), na.rm = TRUE), .groups = "drop")

  ## Add rt flag----
  dat_flag <- dat_flag %>%
    left_join(dat_rt_stds, by = "precursor_ion_name") %>%
    mutate(rt_flag = ifelse((abs(
      retention_time - rt_std
    ) > rt_flex),
    "rt_flag", NA
    ))

  ## Get area of blanks-----
  dat_area_blanks <- dat_og %>%
    filter(str_detect(replicate_name, blank_flag)) %>%
    select(precursor_ion_name, area) %>%
    group_by(precursor_ion_name) %>%
    summarise(area_blank = mean((area), na.rm = TRUE), .groups = "drop")

  ## Add area blank flag----
  dat_flag <- dat_flag %>%
    left_join(dat_area_blanks, by = "precursor_ion_name") %>%
    mutate(blank_flag = ifelse(area / area_blank < blank_ratio_max, "blank_flag", NA))

  # Finally, combine all the flags and throw out any peak with a flag
  dat_flag <- dat_flag %>%
    mutate(flags = paste(sn_flag, ppm_flag, areamin_flag, rt_flag, blank_flag,
      sep = ", "
    )) %>%
    mutate(flags = as.character(flags %>%
      str_remove_all("NA, ") %>%
      str_remove_all("NA"))) %>%
    mutate(QC_area = ifelse(str_detect(flags, "flag"), NA, area))

  # Combine with cleaned up data ----
  dat_final <- dat_og %>%
    left_join(
      dat_flag %>%
        select(precursor_ion_name, replicate_name, flags, QC_area),
      by = c("replicate_name", "precursor_ion_name")
    )

  # To write the file
  comment_text <- paste0(
    "Hello! Here are the QC parameters ",
    "Minimum area for a real peak: ",
    area_min,
    ". ",
    "RT flexibility: ",
    rt_flex,
    " minutes. Blank can be this fraction of a sample: ",
    blank_ratio_max,
    ". S/N ratio: ",
    sn_min,
    ". Parts per million flexibility: ",
    ppm_flex,
    ". Processed on: ",
    Sys.time(),
    ". "
  )
  new_filename <- fileout
  con <- file(new_filename, open = "wt")
  writeLines(paste(comment_text), con)
  write.csv(dat_final, con)
  close(con)
  print(paste("Success! File written to ", fileout))
}

# BMIS function ------
BMIS <-
  function(sample_key,
           is_names_file,
           dat_pos_file,
           dat_neg_file = NA,
           cut_off1,
           cut_off2,
           smps_to_dump = c(),
           is_to_dump = c()) {
    # Import data -----
    sample_key_all <- read_csv(sample_key, show_col_types = FALSE)
    is_names <- read_csv(is_names_file, show_col_types = FALSE)

    # Tidy up samp and pooled data -----
    dat_pos <-
      read_csv(dat_pos_file, skip = 1, show_col_types = FALSE) %>%
      mutate(fraction = "HILICPos")
    if (!is.na(dat_neg_file)) {
      dat_neg <-
        read_csv(dat_neg_file, skip = 1, show_col_types = FALSE) %>%
        mutate(fraction = "HILICNeg")
      dat <- bind_rows(dat_pos, dat_neg)
    } else {
      dat <- dat_pos
    }

    # drop blanks, standards, bad IS, and bad samps
    dat <- dat %>%
      filter(!str_detect(replicate_name, "_Blk_")) %>%
      filter(!str_detect(replicate_name, "_Std_")) %>%
      mutate(area = as.numeric(area)) %>%
      filter(!(replicate_name %in% smps_to_dump)) %>%
      filter(!(precursor_ion_name %in% is_to_dump))


    # Tidy up internal standard data ----
    is_dat_full <- dat %>%
      filter(precursor_ion_name %in% is_names$is_name)
    is_dat <- is_dat_full %>%
      select(replicate_name, precursor_ion_name, area) %>%
      rename(mass_feature = precursor_ion_name)

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

    # Remove is data from dat file
    dat <- dat %>%
      filter(!precursor_ion_name %in% is_names$is_name)

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
      rename(mass_feature = precursor_ion_name) %>%
      select(replicate_name, mass_feature, fraction, area, QC_area) %>%
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
    for (i in 1:length(unique(is_dat$mass_feature))) {
      split_dat[[i]] <-
        # is_dat %>%
        # mutate(QC_area = area) %>%
        bind_rows(is_dat %>%
          mutate(QC_area = area), dat) %>%
        mutate(MIS = unique(is_dat$mass_feature)[i]) %>%
        left_join(
          is_dat %>%
            rename(MIS = mass_feature, is_area = area) %>%
            select(MIS, replicate_name, is_area),
          by = c("replicate_name", "MIS")
        ) %>%
        left_join(
          is_means %>%
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

    poodat <- poodat %>% left_join(
      poodat %>%
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
      round(length(Try$mass_feature) / length(newpoodat$mass_feature), digits = 3) * 100,
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

# Targeted search function -------
ms1_search <-
  function(search_compound_filename,
           folder_to_search,
           out_file_name,
           ion_mode = "positive",
           sample_flag,
           blank_flag,
           meta_data_file,
           ppm_flex = 5) {
    ## get config data ----
    cmps_to_search <- read_csv(
      search_compound_filename,
      show_col_types = FALSE
    )
    files <- list.files(folder_to_search, pattern = ".mzML") %>%
      str_subset(paste0("(", sample_flag, ")|(", blank_flag, ")"))
    meta_dat <- read_csv(
      meta_data_file,
      show_col_types = FALSE
    ) %>%
      mutate(
        filename = paste0(replicate_name, ".mzML")
      )

    ## get ms1 data ----
    print("Grabbing MS1 data from files")
    ms1data <- grabMSdata(
      files = paste0(folder_to_search, "/", files),
      grab_what = c("MS1")
    )

    # prep figure ------
    pdf(out_file_name,
      width = 11,
      height = 8
    )

    # make each compound plot ----
    for (i in 1:length(cmps_to_search$compound_name)) {
      compound_name <- cmps_to_search$compound_name[i]
      if (ion_mode == "negative") {
        mz_oi <- cmps_to_search$mz_neg[i]
      } else if (ion_mode == "positive") {
        mz_oi <- cmps_to_search$mz_pos[i]
      }

      print(paste("Plotting", compound_name))

      rts <- ms1data$MS1 %>%
        select(rt, filename) %>%
        distinct()

      ms1dat_toplot <- ms1data$MS1 %>%
        filter(mz %between% pmppm(mz_oi, ppm_flex)) %>%
        group_by(rt, filename) %>%
        summarise(int = sum(int)) %>%
        ungroup() %>%
        right_join(rts, by = c("rt", "filename")) %>%
        mutate(
          int = ifelse(is.na(int), 0, int),
          target_mass = mz_oi
        ) %>%
        ungroup() %>%
        left_join(
          meta_dat %>%
            select(filename, treatment),
          by = "filename"
        ) %>%
        mutate(treatment = ifelse(str_detect(filename, blank_flag), "Blank", treatment))

      ms1dat_toplot <- ms1dat_toplot %>%
        arrange(rt, filename) %>%
        group_by(filename) %>%
        # if median of 3 = 0, only one had a signal, replace with 0
        mutate(int_roll = rollapply(int, 3, median, align = "right", fill = 0)) %>%
        mutate(int = ifelse(int_roll == 0, 0, int)) %>%
        ungroup()


      ## construct plot -----
      g <- ggplot(
        data = ms1dat_toplot,
        aes(
          x = rt,
          y = int_roll,
          group = filename
        )
      ) +
        geom_line() +
        facet_wrap(
          facets = vars(treatment),
          ncol = 1,
          # scale = "free_y",
          strip.position = "right"
        ) +
        labs(
          title = paste0(compound_name, ", ppm = ", ppm_flex),
          subtitle = paste0("mass = ", round(mz_oi, digits = 4)),
          y = "Intensity",
          x = "Retention time, (min)"
        ) +
        theme_bw() +
        theme(
          strip.background = element_blank(),
          panel.spacing = unit(0, "lines")
        )
      plot(g) %>% suppressWarnings()
    }

    dev.off()
  }
