calc_etr_point <- function() {

}

calc_etr_spatial <- function() {
  list.files(
    "~/MCO_onedrive/General/nexgddp_cmip6_montana/data-derived/nexgddp_cmip6/",
    full.names = T,
    pattern = "ssp585"
  ) %>%
    grep(".aux.json", ., value = T, invert = T) %>%
    grep("MRI-ESM2-0", ., value = T)
}
