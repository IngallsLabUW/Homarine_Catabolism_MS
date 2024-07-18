#' Make EIC plot for mz, rt pair of interest
#'
#' @param msdata dataframe of ms1 data with the following columns: mz, rt, int, filename
#' @param mz mz value of interest
#' @param rt rt value of interest (in minutes)
#' @param samples_oi dataframe of samples of interest with the following columns: filename, treatment
#' @param mz_ppm mz ppm tolerance. Default is 5
#' @param rt_buffer rt buffer in minutes. Default is 2.5
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_EIC <- function(ms1data, m_z, r_t, samples_oi, mz_ppm = 5, rt_buffer = 2.5){
    rt_dat <- ms1data %>% select(rt, filename) %>% distinct() %>%
        filter(rt %between% c(r_t-rt_buffer, r_t+rt_buffer))

    ms1_data <- ms1data[mz %between% pmppm(m_z, mz_ppm) & rt %between%c(r_t-rt_buffer, r_t+rt_buffer)] %>%
        full_join(rt_dat, by = join_by(rt, filename)) %>%
        left_join(samples_oi %>%
                      select(filename, treatment),
                  by = join_by(filename)) %>%
        mutate(int = ifelse(is.na(int), 0, int))

    suppressWarnings(g_ms1 <- ggplot(ms1_data) +
                             geom_line(aes(x = rt, y = int, group  = filename)) +
                             facet_wrap(~treatment, ncol = 3) +
                             theme_bw() +
                             theme(legend.position = "none") +
                             geom_vline("vline", xintercept = r_t, linetype = "dashed") +
                             scale_y_continuous(expand = c(0, NA)) +
                             labs(y = "Intensity",
                                  x = "Retention Time (min)") )
    return(g_ms1)
}

plot_EIC2 <- function(ms1data, m_z, r_t, samples_oi, mz_ppm = 5, rt_buffer = 2.5, scaler = NA){
    rt_dat <- ms1data %>% select(rt, filename) %>% distinct() %>%
        filter(rt %between% c(r_t-rt_buffer, r_t+rt_buffer))

    ms1_data <- ms1data[mz %between% pmppm(m_z, mz_ppm) & rt %between%c(r_t-rt_buffer, r_t+rt_buffer)] %>%
        full_join(rt_dat, by = join_by(rt, filename)) %>%
        left_join(samples_oi %>%
                      select(filename, treatment_short, experiment),
                  by = join_by(filename)) %>%
        mutate(int = ifelse(is.na(int), 0, int))
    if (!is.na(scaler)){
        ms1_data <- ms1_data %>% mutate(int = int*scaler)
    }
    suppressWarnings(g_ms1 <- ggplot(ms1_data) +
                         geom_line(aes(x = rt, y = int, group  = filename, color = treatment_short)) +
                         facet_wrap(~experiment, scales = "free_y") +
                         theme_bw() +
                         #theme(legend.position = "none") +
                         #geom_vline("vline", xintercept = r_t, linetype = "dashed") +
                         scale_y_continuous(expand = c(0, NA)) +
                         labs(y = "Intensity",
                              x = "Retention Time (min)") )
    return(g_ms1)
}

#' Plot spectrum
#'
#' @param spec_data dataframe of spectrum data with the following columns: mz_round, int_norm.
#' @return ggplot object
#'
#' @note This function is used to plot the spectrum of a given mz, rt pair after the data have been normalized (to max of 1) and mzs are binned
plot_spectrum <- function(spec_data){
    g_ms <- ggplot(ms2data) +
        geom_segment(
            aes(
                x = mz_round,
                xend = mz_round,
                y = 0,
                yend = int_norm
            ),
            alpha = 0.5,
        ) +
        geom_text_repel(
            aes(
                x = mz_round,
                y = int_norm,
                label = mz_round
            ),
            nudge_y = 0.05,
            size = 2.5
        ) +
        theme_bw()+
        labs(
            y = NULL,
            x = "m/z") +
        scale_y_continuous(expand = c(0, NA),
                           limits = c(0, 1.07)) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    return(g_ms)

}

#' Pull MS2 data for a given mz, rt pair
#'
#' @param dda_data dataframe of ms2 data with the following columns: premz, fragmz, rt, int, filename
#' @param m_z mz value of interest (precursor mz)
#' @param r_t rt value of interest (in minutes)
#' @param mz_ppm mz ppm tolerance. Default is 5
#' @param rt_buffer rt buffer in minutes. Default is 0.2
#'
pull_ms2_data <- function(dda_data, m_z, r_t, mz_ppm = 5, rt_buffer = 0.2, mz_round = 3){
ms2data <-  dda_data[premz %between% pmppm(m_z, mz_ppm) & rt %between%c(r_t-rt_buffer, r_t+rt_buffer)] %>%
    mutate(mz_round = round(fragmz, digits = 3)) %>%
    group_by(mz_round) %>%
    summarise(int = sum(int)) %>%
    ungroup() %>%
    mutate(int_norm = int/max(int)) %>%
    filter(int_norm > 0.01)
return(ms2data)

# ms2data <-  msdata_dda$MS2[premz %between% pmppm(mf_info$mz, 50) & rt %between%c(mf_info$RT-1, mf_info$RT+1)] %>%
#     select(fragmz, int) %>%
#     arrange(fragmz) %>%
#     mutate(lag_flag = (fragmz - lag(fragmz)>0.001) | (fragmz - lag(fragmz, -1)>0.001))

}