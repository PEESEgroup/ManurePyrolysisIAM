# Copyright 2019 Battelle Memorial Institute; see the LICENSE file.

#' module_aglu_L142.ag_Fert_IO_R_C_Y_GLU
#'
#' Calculate the adjusted fertilizer production by country / year, fertilizer net exports by GCAM region / year,
#' and fertilizer input-output coefficients by GCAM region / commodity / year / GLU.
#'
#' @param command API command to execute
#' @param ... other optional parameters, depending on command
#' @return Depends on \code{command}: either a vector of required inputs,
#' a vector of output names, or (if \code{command} is "MAKE") all
#' the generated outputs: \code{L142.ag_Fert_Prod_MtN_ctry_Y}, \code{L142.ag_Fert_NetExp_MtN_R_Y}, \code{L142.ag_Fert_IO_R_C_Y_GLU}. The corresponding file in the
#' original data system was \code{LB142.ag_Fert_IO_R_C_Y_GLU.R} (aglu level1).
#' @details This chunk calculates fertilizer production by country / year (adjusted to global total consumption),
#' fertilizer net exports by GCAM region / year as production minus consumption, and fertilizer input-output coefficients
#' by GCAM region / commodity / year / GLU.
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter full_join group_by left_join mutate right_join select semi_join summarise
#' @importFrom tidyr complete replace_na
#' @author RC June 2017
module_aglu_L142.ag_Fert_IO_R_C_Y_GLU <- function(command, ...) {

  MODULE_INPUTS <-
    c(FILE = "common/iso_GCAM_regID",
      FILE = "common/GCAM_region_names",
      FILE = "aglu/FAO/FAO_ag_items_PRODSTAT",
      FILE = "aglu/A_SSP1_RCP2p6_an_prod_supply",
      "L100.LDS_ag_prod_t",
      "L100.FAO_Fert_Cons_tN",
      "L100.FAO_Fert_Prod_tN",
      "L100.FAO_Fert_Cons_tK2O",
      "L100.FAO_Fert_Prod_tK2O",
      "L100.FAO_Fert_Cons_tP2O5",
      "L100.FAO_Fert_Prod_tP2O5",
      "L101.ag_Prod_Mt_R_C_Y_GLU",
      "L141.ag_Fert_Cons_MtN_ctry_crop",
      "L141.ag_Fert_Cons_MtK2O_ctry_crop",
      "L141.ag_Fert_Cons_MtP2O5_ctry_crop")

  MODULE_OUTPUTS <-
    c("L142.ag_Fert_Prod_MtN_ctry_Y",
      "L142.ag_Fert_NetExp_MtN_R_Y",
      "L142.ag_Fert_IO_R_C_Y_GLU",
      "L142.ag_Fert_IO_R_C_Y_GLU_K2O",
      "L142.ag_Fert_IO_R_C_Y_GLU_P2O5",
      "L142.ag_Fert_IO_R_C_Y_GLU_biochar")

  if(command == driver.DECLARE_INPUTS) {
    return(MODULE_INPUTS)
  } else if(command == driver.DECLARE_OUTPUTS) {
    return(MODULE_OUTPUTS)
  } else if(command == driver.MAKE) {

    Fert_Cons_MtN <- Fert_Cons_MtN_unscaled <- Fert_IO <- Fert_IO_unscaled <- Prod_share <-
      prod <- cons <- total <- adj <- scaler <- GCAM_commodity <- GCAM_region_ID <- GTAP_crop <-
      GLU <- iso <- value <- year <- GCAM_subsector <- NULL   # silence package checks

    all_data <- list(...)[[1]]

    # Load required inputs ----
    get_data_list(all_data, MODULE_INPUTS, strip_attributes = TRUE)


    # Compile N fertilizer production and consumption by country, and adjust country production so that production and consumption balance globally
    L100.FAO_Fert_Prod_tN %>%
      select(iso, year, prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tN, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      group_by(year) %>%
      # Calculate the global total production and consumption
      summarise(prod = sum(prod), cons = sum(cons)) %>%
      ungroup() %>%
      # Calculate the rate to adjust production so that global production equals consumption
      mutate(adj = cons / prod) %>%
      select(year, adj) ->
      L142.ag_Fert_Prod_adj_N

    L100.FAO_Fert_Prod_tK2O %>%
      select(iso, year, prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tK2O, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      group_by(year) %>%
      # Calculate the global total production and consumption
      summarise(prod = sum(prod), cons = sum(cons)) %>%
      ungroup() %>%
      # Calculate the rate to adjust production so that global production equals consumption
      mutate(adj = cons / prod) %>%
      select(year, adj) ->
      L142.ag_Fert_Prod_adj_K2O

    L100.FAO_Fert_Prod_tP2O5 %>%
      select(iso, year, prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tP2O5, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      group_by(year) %>%
      # Calculate the global total production and consumption
      summarise(prod = sum(prod), cons = sum(cons)) %>%
      ungroup() %>%
      # Calculate the rate to adjust production so that global production equals consumption
      mutate(adj = cons / prod) %>%
      select(year, adj) ->
      L142.ag_Fert_Prod_adj_P2O5

    L100.FAO_Fert_Prod_tN %>%
      select(iso, year, value) %>%
      left_join_error_no_match(L142.ag_Fert_Prod_adj_N, by = "year") %>%   # Match in the rates for adjustment
      mutate(value = value * adj,                                        # Adjust production
             value = value * CONV_T_MT) %>%                              # Convert unit of production from tons to million tons of Nitrogen
      select(-adj) ->
      L142.ag_Fert_Prod_MtN_ctry_Y

    L100.FAO_Fert_Prod_tK2O %>%
      select(iso, year, value) %>%
      left_join_error_no_match(L142.ag_Fert_Prod_adj_K2O, by = "year") %>%   # Match in the rates for adjustment
      mutate(value = value * adj,                                        # Adjust production
             value = value * CONV_T_MT) %>%                              # Convert unit of production from tons to million tons of Nitrogen
      select(-adj) ->
      L142.ag_Fert_Prod_MtK2O_ctry_Y

    L100.FAO_Fert_Prod_tP2O5 %>%
      select(iso, year, value) %>%
      left_join_error_no_match(L142.ag_Fert_Prod_adj_P2O5, by = "year") %>%   # Match in the rates for adjustment
      mutate(value = value * adj,                                        # Adjust production
             value = value * CONV_T_MT) %>%                              # Convert unit of production from tons to million tons of Nitrogen
      select(-adj) ->
      L142.ag_Fert_Prod_MtP2O5_ctry_Y

    # Aggregate N fertilizer adjusted production and consumption to GCAM region level to calculate net exports
    L142.ag_Fert_Prod_MtN_ctry_Y %>%
      rename(prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tN, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%           # Match in GCAM region ID
      group_by(GCAM_region_ID, year) %>%
      summarise(prod = sum(prod), cons = sum(cons)) %>%                  # Aggregate to region total
      ungroup() %>%                                                      # Ungroup before complete
      mutate(cons = cons * CONV_T_MT,                                    # Convert unit of consumption from tons to million tons of Nitrogen
             value = prod - cons,                                        # Calculate net exports as production minus consumption
             GCAM_commodity = aglu.FERT_NAME) %>%                        # Add GCAM commodity category for N fertilizer
      select(-prod) %>%                                                  # Only regional consumption and net exports are needed
      complete(GCAM_region_ID = unique(iso_GCAM_regID$GCAM_region_ID),
               GCAM_commodity, year, fill = list(cons = 0, value = 0)) ->  # Fill in missing regions with 0
      L142.ag_Fert_MtN_R_Y

    L142.ag_Fert_Prod_MtK2O_ctry_Y %>%
      rename(prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tK2O, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%           # Match in GCAM region ID
      group_by(GCAM_region_ID, year) %>%
      summarise(prod = sum(prod), cons = sum(cons)) %>%                  # Aggregate to region total
      ungroup() %>%                                                      # Ungroup before complete
      mutate(cons = cons * CONV_T_MT,                                    # Convert unit of consumption from tons to million tons of Nitrogen
             value = prod - cons,                                        # Calculate net exports as production minus consumption
             GCAM_commodity = aglu.FERT_NAME) %>%                        # Add GCAM commodity category for N fertilizer
      select(-prod) %>%                                                  # Only regional consumption and net exports are needed
      complete(GCAM_region_ID = unique(iso_GCAM_regID$GCAM_region_ID),
               GCAM_commodity, year, fill = list(cons = 0, value = 0)) ->  # Fill in missing regions with 0
      L142.ag_Fert_MtK2O_R_Y

    L142.ag_Fert_Prod_MtP2O5_ctry_Y %>%
      rename(prod = value) %>%
      # Combine with fertilizer consumption, use full_join to keep all observations, such as ones only have consumption
      full_join(select(L100.FAO_Fert_Cons_tP2O5, iso, year, cons = value), by = c("iso", "year")) %>%
      replace_na(list(prod = 0, cons = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%           # Match in GCAM region ID
      group_by(GCAM_region_ID, year) %>%
      summarise(prod = sum(prod), cons = sum(cons)) %>%                  # Aggregate to region total
      ungroup() %>%                                                      # Ungroup before complete
      mutate(cons = cons * CONV_T_MT,                                    # Convert unit of consumption from tons to million tons of Nitrogen
             value = prod - cons,                                        # Calculate net exports as production minus consumption
             GCAM_commodity = aglu.FERT_NAME) %>%                        # Add GCAM commodity category for N fertilizer
      select(-prod) %>%                                                  # Only regional consumption and net exports are needed
      complete(GCAM_region_ID = unique(iso_GCAM_regID$GCAM_region_ID),
               GCAM_commodity, year, fill = list(cons = 0, value = 0)) ->  # Fill in missing regions with 0
      L142.ag_Fert_MtP2O5_R_Y

    # Separate the table for consumption by region / year
    L142.ag_Fert_MtN_R_Y %>%
      select(-value) ->
      L142.ag_Fert_Cons_MtN_R_Y

    L142.ag_Fert_MtK2O_R_Y %>%
      select(-value) ->
      L142.ag_Fert_Cons_MtK2O_R_Y

    L142.ag_Fert_MtP2O5_R_Y %>%
      select(-value) ->
      L142.ag_Fert_Cons_MtP2O5_R_Y

    # Separate the table for net exports by region / year
    L142.ag_Fert_MtN_R_Y %>%
      select(-cons) ->
      L142.ag_Fert_NetExp_MtN_R_Y

    # Calculate fertilizer input-output coefficients, scaled so that consumption of fertilizer balance
    # First, downscale fertilizer demands by country and crop to GLU
    # NOTE: Allocate fertilizer consumption to GLUs on the basis of production, not harvested area
    # Calculate agriculture prodcution total by country and crop
    L100.LDS_ag_prod_t %>%
      group_by(iso, GTAP_crop) %>%
      summarise(total = sum(value)) %>%
      ungroup() ->
      L142.ag_Prod_t_ctry_crop

    # Start with agricultural production by country / crop / GLU, calculate the production share of GLU to country,
    # Use the share to downscale country fertilizer demand to GLU.
    L100.LDS_ag_prod_t %>%
      # Match in country total production
      left_join_error_no_match(L142.ag_Prod_t_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Match in country fertilizer consumption by crop
      left_join(L141.ag_Fert_Cons_MtN_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Calculate production share of GLU to country total
      mutate(Prod_share = value / total,
             # Downscale fertilizer demands to GLU using the production share
             Fert_Cons_MtN = Fert_Cons_MtN * Prod_share) %>%
      replace_na(list(Fert_Cons_MtN = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%
      # Map in GCAM commodity, creates NA, use left_join instead of left_join_error_no_match
      left_join(select(FAO_ag_items_PRODSTAT, GTAP_crop, GCAM_commodity, GCAM_subsector), by = "GTAP_crop") %>%
      # Drop crops not belong to GCAM commodity
      filter(!is.na(GCAM_commodity)) %>%
      # Aggregate fertilizer demands by GCAM region, commodity, and GLU
      group_by(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU) %>%
      summarise(Fert_Cons_MtN = sum(Fert_Cons_MtN)) %>%
      ungroup() %>%
      # Match in agricultural production by GCAM region / commodity / GLU in the base year; this creates NAs
      left_join(filter(L101.ag_Prod_Mt_R_C_Y_GLU, year %in% aglu.BASE_YEAR_IFA),
                by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled input-output coefficients as unscaled fertilizer demands divided by agricultural production
      mutate(Fert_IO_unscaled = Fert_Cons_MtN / value,
             Fert_IO_unscaled = replace(Fert_IO_unscaled, Fert_IO_unscaled == Inf, 0)) %>%
      select(-year, -value, -Fert_Cons_MtN) %>%
      # Match these coefficients into historical agricultural production (right_join: same coefficients for all years)
      right_join(L101.ag_Prod_Mt_R_C_Y_GLU, by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled fertilizer consumption by year as production multiplied by input-output coefficients
      # GPK note 09/2018 - moving the replace_na() command from a few lines up to here, in order to accommodate missing
      # values which can occur from country/crop observations in FAOSTAT but not Monfreda/LDS (e.g., Puerto Rico rice)
      replace_na(list(Fert_IO_unscaled = 0)) %>%
      mutate(Fert_Cons_MtN_unscaled = value * Fert_IO_unscaled) ->
      L142.ag_Fert_Cons_MtN_R_C_Y_GLU

    # Compute region/year scalers so that sum of fertilizer consumption for all commodities equals total consumption
    L142.ag_Fert_Cons_MtN_R_C_Y_GLU %>%
      group_by(GCAM_region_ID, year) %>%
      # Caclulate total unscaled regional consumption
      summarise(Fert_Cons_MtN_unscaled = sum(Fert_Cons_MtN_unscaled)) %>%
      ungroup() %>%
      # Match in historical total consumption by region
      left_join_error_no_match(L142.ag_Fert_Cons_MtN_R_Y, by = c("GCAM_region_ID", "year")) %>%
      # Calculate the regional scaler so that consumption balance
      mutate(scaler = cons / Fert_Cons_MtN_unscaled) %>%
      replace_na(list(scaler = 0)) %>%
      select(GCAM_region_ID, year, scaler) %>%
      # Match the scalers to unscaled fertilizer consumption by region/commodity/GLU (right_join: same scaler for individual region)
      right_join(L142.ag_Fert_Cons_MtN_R_C_Y_GLU, by = c("GCAM_region_ID", "year")) %>%
      # Calculate scaled consumption
      mutate(Fert_Cons_MtN = Fert_Cons_MtN_unscaled * scaler,
             # Calculate the scalced input-output coefficient
             Fert_IO = Fert_IO_unscaled * scaler) %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, value = Fert_IO) ->
      L142.ag_Fert_IO_R_C_Y_GLU



    # Start with agricultural production by country / crop / GLU, calculate the production share of GLU to country,
    # Use the share to downscale country fertilizer demand to GLU.
    L100.LDS_ag_prod_t %>%
      # Match in country total production
      left_join_error_no_match(L142.ag_Prod_t_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Match in country fertilizer consumption by crop
      left_join(L141.ag_Fert_Cons_MtK2O_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Calculate production share of GLU to country total
      mutate(Prod_share = value / total,
             # Downscale fertilizer demands to GLU using the production share
             Fert_Cons_MtK2O = Fert_Cons_MtK2O * Prod_share) %>%
      replace_na(list(Fert_Cons_MtK2O = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%
      # Map in GCAM commodity, creates NA, use left_join instead of left_join_error_no_match
      left_join(select(FAO_ag_items_PRODSTAT, GTAP_crop, GCAM_commodity, GCAM_subsector), by = "GTAP_crop") %>%
      # Drop crops not belong to GCAM commodity
      filter(!is.na(GCAM_commodity)) %>%
      # Aggregate fertilizer demands by GCAM region, commodity, and GLU
      group_by(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU) %>%
      summarise(Fert_Cons_MtK2O = sum(Fert_Cons_MtK2O)) %>%
      ungroup() %>%
      # Match in agricultural production by GCAM region / commodity / GLU in the base year; this creates NAs
      left_join(filter(L101.ag_Prod_Mt_R_C_Y_GLU, year %in% aglu.BASE_YEAR_IFA),
                by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled input-output coefficients as unscaled fertilizer demands divided by agricultural production
      mutate(Fert_IO_unscaled = Fert_Cons_MtK2O / value,
             Fert_IO_unscaled = replace(Fert_IO_unscaled, Fert_IO_unscaled == Inf, 0)) %>%
      select(-year, -value, -Fert_Cons_MtK2O) %>%
      # Match these coefficients into historical agricultural production (right_join: same coefficients for all years)
      right_join(L101.ag_Prod_Mt_R_C_Y_GLU, by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled fertilizer consumption by year as production multiplied by input-output coefficients
      # GPK note 09/2018 - moving the replace_na() command from a few lines up to here, in order to accommodate missing
      # values which can occur from country/crop observations in FAOSTAT but not Monfreda/LDS (e.g., Puerto Rico rice)
      replace_na(list(Fert_IO_unscaled = 0)) %>%
      mutate(Fert_Cons_MtK2O_unscaled = value * Fert_IO_unscaled) ->
      L142.ag_Fert_Cons_MtK2O_R_C_Y_GLU

    # Compute region/year scalers so that sum of fertilizer consumption for all commodities equals total consumption
    L142.ag_Fert_Cons_MtK2O_R_C_Y_GLU %>%
      group_by(GCAM_region_ID, year) %>%
      # Caclulate total unscaled regional consumption
      summarise(Fert_Cons_MtK2O_unscaled = sum(Fert_Cons_MtK2O_unscaled)) %>%
      ungroup() %>%
      # Match in historical total consumption by region
      left_join_error_no_match(L142.ag_Fert_Cons_MtK2O_R_Y, by = c("GCAM_region_ID", "year")) %>%
      # Calculate the regional scaler so that consumption balance
      mutate(scaler = cons / Fert_Cons_MtK2O_unscaled) %>%
      replace_na(list(scaler = 0)) %>%
      select(GCAM_region_ID, year, scaler) %>%
      # Match the scalers to unscaled fertilizer consumption by region/commodity/GLU (right_join: same scaler for individual region)
      right_join(L142.ag_Fert_Cons_MtK2O_R_C_Y_GLU, by = c("GCAM_region_ID", "year")) %>%
      # Calculate scaled consumption
      mutate(Fert_Cons_MtK2O = Fert_Cons_MtK2O_unscaled * scaler,
             # Calculate the scalced input-output coefficient
             Fert_IO = Fert_IO_unscaled * scaler) %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, value = Fert_IO) ->
      L142.ag_Fert_IO_R_C_Y_GLU_K2O


    # Start with agricultural production by country / crop / GLU, calculate the production share of GLU to country,
    # Use the share to downscale country fertilizer demand to GLU.
    L100.LDS_ag_prod_t %>%
      # Match in country total production
      left_join_error_no_match(L142.ag_Prod_t_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Match in country fertilizer consumption by crop
      left_join(L141.ag_Fert_Cons_MtP2O5_ctry_crop, by = c("iso", "GTAP_crop")) %>%
      # Calculate production share of GLU to country total
      mutate(Prod_share = value / total,
             # Downscale fertilizer demands to GLU using the production share
             Fert_Cons_MtP2O5 = Fert_Cons_MtP2O5 * Prod_share) %>%
      replace_na(list(Fert_Cons_MtN = 0)) %>%
      left_join_error_no_match(iso_GCAM_regID, by = "iso") %>%
      # Map in GCAM commodity, creates NA, use left_join instead of left_join_error_no_match
      left_join(select(FAO_ag_items_PRODSTAT, GTAP_crop, GCAM_commodity, GCAM_subsector), by = "GTAP_crop") %>%
      # Drop crops not belong to GCAM commodity
      filter(!is.na(GCAM_commodity)) %>%
      # Aggregate fertilizer demands by GCAM region, commodity, and GLU
      group_by(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU) %>%
      summarise(Fert_Cons_MtP2O5 = sum(Fert_Cons_MtP2O5)) %>%
      ungroup() %>%
      # Match in agricultural production by GCAM region / commodity / GLU in the base year; this creates NAs
      left_join(filter(L101.ag_Prod_Mt_R_C_Y_GLU, year %in% aglu.BASE_YEAR_IFA),
                by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled input-output coefficients as unscaled fertilizer demands divided by agricultural production
      mutate(Fert_IO_unscaled = Fert_Cons_MtP2O5 / value,
             Fert_IO_unscaled = replace(Fert_IO_unscaled, Fert_IO_unscaled == Inf, 0)) %>%
      select(-year, -value, -Fert_Cons_MtP2O5) %>%
      # Match these coefficients into historical agricultural production (right_join: same coefficients for all years)
      right_join(L101.ag_Prod_Mt_R_C_Y_GLU, by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) %>%
      # Calculate unscaled fertilizer consumption by year as production multiplied by input-output coefficients
      # GPK note 09/2018 - moving the replace_na() command from a few lines up to here, in order to accommodate missing
      # values which can occur from country/crop observations in FAOSTAT but not Monfreda/LDS (e.g., Puerto Rico rice)
      replace_na(list(Fert_IO_unscaled = 0)) %>%
      mutate(Fert_Cons_MtP2O5_unscaled = value * Fert_IO_unscaled) ->
      L142.ag_Fert_Cons_MtP2O5_R_C_Y_GLU

    # Compute region/year scalers so that sum of fertilizer consumption for all commodities equals total consumption
    L142.ag_Fert_Cons_MtP2O5_R_C_Y_GLU %>%
      group_by(GCAM_region_ID, year) %>%
      # Caclulate total unscaled regional consumption
      summarise(Fert_Cons_MtP2O5_unscaled = sum(Fert_Cons_MtP2O5_unscaled)) %>%
      ungroup() %>%
      # Match in historical total consumption by region
      left_join_error_no_match(L142.ag_Fert_Cons_MtP2O5_R_Y, by = c("GCAM_region_ID", "year")) %>%
      # Calculate the regional scaler so that consumption balance
      mutate(scaler = cons / Fert_Cons_MtP2O5_unscaled) %>%
      replace_na(list(scaler = 0)) %>%
      select(GCAM_region_ID, year, scaler) %>%
      # Match the scalers to unscaled fertilizer consumption by region/commodity/GLU (right_join: same scaler for individual region)
      right_join(L142.ag_Fert_Cons_MtP2O5_R_C_Y_GLU, by = c("GCAM_region_ID", "year")) %>%
      # Calculate scaled consumption
      mutate(Fert_Cons_MtP2O5 = Fert_Cons_MtP2O5_unscaled * scaler,
             # Calculate the scalced input-output coefficient
             Fert_IO = Fert_IO_unscaled * scaler) %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, value = Fert_IO) ->
      L142.ag_Fert_IO_R_C_Y_GLU_P2O5


    # Check to make sure that the fertilizer inputs do not blink in and out (if present in any year, need to be present in all years)
    L142.ag_Fert_IO_R_C_Y_GLU %>%
      group_by(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU) %>%
      summarise(value = sum(value)) %>%                 # Get the total of all years
      ungroup() %>%
      filter(value != 0) %>%                            # Filter the region/commodity/GLU that are not completely missing for all years
      select(-value) %>%
      unique() ->                                       # Select the region/commodity/GLU combinations presented
      L142.Fert_IO_check
    L142.ag_Fert_IO_R_C_Y_GLU %>%
      # Filter the observations with the selected region/commodity/GLU combinations
      semi_join(L142.Fert_IO_check, by = c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU")) ->
      L142.Fert_IO_check

    # For those region/commodity/GLU that are not completely missing for all years, no missing for any year
    if(any(L142.Fert_IO_check$value == 0)) {
      stop("Fertilizer input-output coefficients need to be specified in all historical years")
    }

    # compared to N, it's fine if K2O and P2O5 are 0, as those types of fertilizers may not be applied

    print(L142.ag_Fert_IO_R_C_Y_GLU_K2O)
    print(L142.ag_Fert_IO_R_C_Y_GLU_P2O5)

    # load in manure production
    # constants for biochar yield
    beef_yield =1/2.105
    dairy_yield =1/2.105
    goat_yield =1/2.055
    pork_yield = 1/2.136
    poultry_yield = 1/2.139

    # constants for P and K nutrients in biochar in kg nutrient/kg biochar
    beef_P = 0.0122
    dairy_P = 0.0122
    goat_P = 0.0053
    pork_P = 0.0750
    poultry_P = 0.0283

    beef_K = 0.0008
    dairy_K = 0.0008
    goat_K = 0.0420
    pork_K = 0.0310
    poultry_K = 0.0746

    K_K2O = 1.2046
    P_P2O5 = 2.2951

    # for each region, calculate the amount of P/K in each representative ton biochar
    print(A_SSP1_RCP2p6_an_prod_supply)
    A_SSP1_RCP2p6_an_prod_supply %>%
      mutate(total_P = if_else(product == "Beef", Yield*beef_yield*beef_P, if_else(product == "Dairy", Yield*dairy_yield*dairy_P, if_else(product == "SheepGoat", Yield*goat_yield*goat_P, if_else(product == "Pork", Yield*pork_yield*pork_P, if_else(product == "Poultry", Yield*poultry_yield*poultry_P, 0))))),
             total_K = if_else(product == "Beef", Yield*beef_yield*beef_K, if_else(product == "Dairy", Yield*dairy_yield*dairy_K, if_else(product == "SheepGoat", Yield*goat_yield*goat_K, if_else(product == "Pork", Yield*pork_yield*pork_K, if_else(product == "Poultry", Yield*poultry_yield*poultry_K, 0))))),
             total_biochar = if_else(product == "Beef", Yield*beef_yield, if_else(product == "Dairy", Yield*dairy_yield, if_else(product == "SheepGoat", Yield*goat_yield, if_else(product == "Pork", Yield*pork_yield, if_else(product == "Poultry", Yield*poultry_yield, 0)))))) %>%
      group_by(GCAM)%>%
      summarise(total_P = sum(total_P), total_K = sum(total_K), total_biochar = sum(total_biochar)) %>%
      ungroup() %>%
      mutate(rep_P2O5 = P_P2O5*total_P/total_biochar, # convert to fertilizer equivalents
             rep_K2O = K_K2O*total_K/total_biochar) %>%
      select(GCAM, rep_P2O5, rep_K2O) -> L142.biochar_nutrient_content # units are in kg K2O/kg biochar

    # match 2015 GCAM data with 2015 fertilizer data
    L142.ag_Fert_IO_R_C_Y_GLU_K2O %>% filter(year==2015) %>%
      left_join_error_no_match(GCAM_region_names, by = "GCAM_region_ID") %>%
      mutate(K2O = value) %>%
      select(-value)-> L142.ag_Fert_IO_R_C_Y_GLU_K2O
    L142.ag_Fert_IO_R_C_Y_GLU_P2O5 %>% filter(year==2015) %>%
      left_join_error_no_match(GCAM_region_names, by = "GCAM_region_ID") %>%
      mutate(P2O5 = value) %>%
      select(-value)->L142.ag_Fert_IO_R_C_Y_GLU_P2O5

    # merge all the data into one table
    print(L142.ag_Fert_IO_R_C_Y_GLU_K2O %>%
            left_join(L142.ag_Fert_IO_R_C_Y_GLU_P2O5, by=c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU", "year", "region")))

    L142.ag_Fert_IO_R_C_Y_GLU_K2O %>%
      left_join(L142.ag_Fert_IO_R_C_Y_GLU_P2O5, by=c("GCAM_region_ID", "GCAM_commodity", "GCAM_subsector", "GLU", "year", "region")) %>%
      left_join(L142.biochar_nutrient_content, by = c("region" = "GCAM")) -> L142_biochar_fertilizer

    # calculate kg biochar per kg crop
    # [kg K2O/ kg crop] / [kg K2O / kg biochar] = kg biochar/kg crop
    L142_biochar_fertilizer %>%
      mutate(kg_biochar_kg_crop_P_limit = K2O/rep_K2O,# above dimensional analysis
             kg_biochar_kg_crop_K_limit = P2O5/rep_P2O5, # above dimensional analysis
             kg_biochar_kg_crop_limited = pmin(kg_biochar_kg_crop_K_limit, kg_biochar_kg_crop_P_limit), # choose whatever limits first
             avoided_K2O = kg_biochar_kg_crop_limited * rep_K2O, # kg biochar/kg crop * kg K2O/kg biochar = kg K2O/kg crop
             avoided_P2O5 = kg_biochar_kg_crop_limited * rep_P2O5) ->
      L142_biochar_fertilizer

    print(L142_biochar_fertilizer)

    L142_biochar_fertilizer %>%
      mutate(value = avoided_K2O) %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, value) ->
      L142.ag_Fert_IO_R_C_Y_GLU_K2O

    L142_biochar_fertilizer %>%
      mutate(value = avoided_P2O5) %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, value) ->
      L142.ag_Fert_IO_R_C_Y_GLU_P2O5

    L142_biochar_fertilizer %>%
      select(GCAM_region_ID, GCAM_commodity, GCAM_subsector, GLU, year, kg_biochar_kg_crop_limited, rep_K2O, rep_P2O5) ->
      L142.ag_Fert_IO_R_C_Y_GLU_biochar


    print(L142.ag_Fert_IO_R_C_Y_GLU_K2O)
    print(L142.ag_Fert_IO_R_C_Y_GLU_P2O5)
    print(L142.ag_Fert_IO_R_C_Y_GLU_biochar)



    # Produce outputs
    L142.ag_Fert_Prod_MtN_ctry_Y %>%
      add_title("Fertilizer production by country / year") %>%
      add_units("Unit = MtN") %>%
      add_comments("Fertilizer production by country is adjusted so that global total production equals consumption") %>%
      add_comments("Units are converted from tons to million tons of Nitrogen") %>%
      add_legacy_name("L142.ag_Fert_Prod_MtN_ctry_Y") %>%
      add_precursors("L100.FAO_Fert_Cons_tN",
                     "L100.FAO_Fert_Prod_tN") ->
      L142.ag_Fert_Prod_MtN_ctry_Y

    L142.ag_Fert_NetExp_MtN_R_Y %>%
      add_title("Fertilizer net exports by GCAM region / year") %>%
      add_units("Unit = MtN") %>%
      add_comments("Fertilizer consumption and adjusted production are aggregated from country to GCAM region level") %>%
      add_comments("Net exports are calculated as production minus consumption, in million tons of Nitrogen") %>%
      add_legacy_name("L142.ag_Fert_NetExp_MtN_R_Y") %>%
      add_precursors("common/iso_GCAM_regID",
                     "L100.FAO_Fert_Cons_tN",
                     "L100.FAO_Fert_Prod_tN") ->
      L142.ag_Fert_NetExp_MtN_R_Y

    L142.ag_Fert_IO_R_C_Y_GLU %>%
      add_title("Fertilizer input-output coefficients by GCAM region / crop / year / GLU") %>%
      add_units("Unitless IO") %>%
      add_comments("Fertilizer demands are downscaled to GLU based on agriculture production share in each country") %>%
      add_comments("Input-output coefficients for each crop are first calculated as fertilizer demands divided by agriculture production in the base year") %>%
      add_comments("And then are scaled so that regional total fertilizer consumptions are balanced") %>%
      add_legacy_name("L142.ag_Fert_IO_R_C_Y_GLU") %>%
      add_precursors("common/iso_GCAM_regID",
                     "aglu/FAO/FAO_ag_items_PRODSTAT",
                     "L100.LDS_ag_prod_t",
                     "L101.ag_Prod_Mt_R_C_Y_GLU",
                     "L141.ag_Fert_Cons_MtN_ctry_crop") ->
      L142.ag_Fert_IO_R_C_Y_GLU

    L142.ag_Fert_IO_R_C_Y_GLU_K2O %>%
      add_title("Fertilizer input-output coefficients by GCAM region / crop / year / GLU") %>%
      add_units("kg avoided K2O/kg crop") %>%
      add_comments("Fertilizer demands are downscaled to GLU based on agriculture production share in each country") %>%
      add_comments("Input-output coefficients for each crop are first calculated as fertilizer demands divided by agriculture production in 2015") %>%
      add_comments("And then are scaled so that regional total fertilizer consumptions are balanced") %>%
      add_legacy_name("L142.ag_Fert_IO_R_C_Y_GLU_K2O") %>%
      add_precursors("common/iso_GCAM_regID",
                     "aglu/FAO/FAO_ag_items_PRODSTAT",
                     "L100.LDS_ag_prod_t",
                     "L101.ag_Prod_Mt_R_C_Y_GLU",
                     "L141.ag_Fert_Cons_MtK2O_ctry_crop") ->
      L142.ag_Fert_IO_R_C_Y_GLU_K2O

    L142.ag_Fert_IO_R_C_Y_GLU_P2O5 %>%
      add_title("Fertilizer input-output coefficients by GCAM region / crop / year / GLU") %>%
      add_units("kg avoided P2O5/kg crop") %>%
      add_comments("Fertilizer demands are downscaled to GLU based on agriculture production share in each country") %>%
      add_comments("Input-output coefficients for each crop are first calculated as fertilizer demands divided by agriculture production in 2015") %>%
      add_comments("And then are scaled so that regional total fertilizer consumptions are balanced") %>%
      add_legacy_name("L142.ag_Fert_IO_R_C_Y_GLU_P2O5") %>%
      add_precursors("common/iso_GCAM_regID",
                     "aglu/FAO/FAO_ag_items_PRODSTAT",
                     "L100.LDS_ag_prod_t",
                     "L101.ag_Prod_Mt_R_C_Y_GLU",
                     "L141.ag_Fert_Cons_MtP2O5_ctry_crop") ->
      L142.ag_Fert_IO_R_C_Y_GLU_P2O5

    L142.ag_Fert_IO_R_C_Y_GLU_biochar %>%
      add_title("Fertilizer input-output coefficients by GCAM region / crop / year / GLU") %>%
      add_units("kg biochar applied/kg crop") %>%
      add_comments("Fertilizer demands are downscaled to GLU based on agriculture production share in each country") %>%
      add_comments("Input-output coefficients for each crop are first calculated as fertilizer demands divided by agriculture production in 2015") %>%
      add_comments("And then are scaled so that regional total fertilizer consumptions are balanced") %>%
      add_legacy_name("L142.ag_Fert_IO_R_C_Y_GLU_biochar") %>%
      add_precursors("common/iso_GCAM_regID",
                     "aglu/FAO/FAO_ag_items_PRODSTAT",
                     "L100.LDS_ag_prod_t",
                     "L101.ag_Prod_Mt_R_C_Y_GLU",
                     "L141.ag_Fert_Cons_MtP2O5_ctry_crop") ->
      L142.ag_Fert_IO_R_C_Y_GLU_biochar

    return_data(MODULE_OUTPUTS)
  } else {
    stop("Unknown command")
  }
}
