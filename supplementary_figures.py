import plotting
import data_manipulation
import constants as c
import pandas as pd


def pop_and_calories(nonBaselineScenario, RCP, SSP, biochar_year):
    """
    plots changes to population and calories consumed
    :param nonBaselineScenario: the scenario to be compared to the released scenario
    :param RCP: the RCP pathways being considered
    :param SSP: the SSP pathways being considered
    :param biochar_year: the year the biochar pathways are being evaluated
    :return: N/A
    """
    # get population data
    released_pop = data_manipulation.get_sensitivity_data(["released"], "population_by_region", SSP, RCP=RCP,
                                                          source="original")
    pyrolysis_pop = data_manipulation.get_sensitivity_data(nonBaselineScenario, "population_by_region", SSP, RCP=RCP,
                                                           source="masked")
    flat_diff_pop = data_manipulation.flat_difference(released_pop, pyrolysis_pop, ["SSP", "GCAM"])
    data_manipulation.drop_missing(flat_diff_pop).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(RCP) + "/flat_change_in_population.csv")

# calculate food accessibility and undernourishment
    released_caloric_consumption = data_manipulation.get_sensitivity_data(["released"], "food_demand_per_capita", SSP,
                                                                          RCP=RCP, source="original")
    pyrolysis_caloric_consumption = data_manipulation.get_sensitivity_data(nonBaselineScenario,
                                                                           "food_demand_per_capita", SSP, RCP=RCP,
                                                                           source="masked")
    released_caloric_consumption = data_manipulation.group(released_caloric_consumption, ["SSP", "GCAM", "Version"])
    pyrolysis_caloric_consumption = data_manipulation.group(pyrolysis_caloric_consumption, ["SSP", "GCAM", "Version"])

    # regional averaged food consumption by food type
    # convert Pcal to kcal/capita/day
    # get population data
    released_pop = data_manipulation.get_sensitivity_data(["released"], "population_by_region", SSP, RCP=RCP,
                                                          source="original")
    pyrolysis_pop = data_manipulation.get_sensitivity_data(nonBaselineScenario, "population_by_region", SSP, RCP=RCP,
                                                           source="masked")
    released_Pcal = data_manipulation.get_sensitivity_data(["released"], "food_consumption_by_type_specific", SSP,
                                                           RCP=RCP, source="original")
    pyrolysis_Pcal = data_manipulation.get_sensitivity_data(nonBaselineScenario, "food_consumption_by_type_specific",
                                                            SSP, RCP=RCP, source="masked")

    # relabel data to make it more human-readable
    released_Pcal['GCAM'] = released_Pcal.apply(lambda row: data_manipulation.relabel_region(row), axis=1)
    released_pop['GCAM'] = released_pop.apply(lambda row: data_manipulation.relabel_region(row), axis=1)
    pyrolysis_Pcal['GCAM'] = pyrolysis_Pcal.apply(lambda row: data_manipulation.relabel_region(row), axis=1)
    pyrolysis_pop['GCAM'] = pyrolysis_pop.apply(lambda row: data_manipulation.relabel_region(row), axis=1)
    released_Pcal['technology'] = released_Pcal.apply(lambda row: data_manipulation.relabel_food(row, "technology"),
                                                      axis=1)
    pyrolysis_Pcal['technology'] = pyrolysis_Pcal.apply(lambda row: data_manipulation.relabel_food(row, "technology"),
                                                        axis=1)

    # drop MiscCrop and FiberCrop because those products don't have meaningful calories and clutter the graph
    released_Pcal = released_Pcal.drop(released_Pcal[released_Pcal["technology"] == "Fiber Crops"].index)
    released_Pcal = released_Pcal.drop(released_Pcal[released_Pcal["technology"] == "Other Crops"].index)
    pyrolysis_Pcal = pyrolysis_Pcal.drop(pyrolysis_Pcal[pyrolysis_Pcal["technology"] == "Fiber Crops"].index)
    pyrolysis_Pcal = pyrolysis_Pcal.drop(pyrolysis_Pcal[pyrolysis_Pcal["technology"] == "Other Crops"].index)

    released_Pcal = data_manipulation.group(released_Pcal, ["GCAM", "SSP", "technology", "Version"])
    pyrolysis_Pcal = data_manipulation.group(pyrolysis_Pcal, ["GCAM", "SSP", "technology", "Version"])

    # calculate food accessibility
    # get food prices
    released_staple_expenditure = data_manipulation.get_sensitivity_data(["released"], "food_demand_prices", SSP,
                                                                         RCP=RCP, source="original")
    pyrolysis_staple_expenditure = data_manipulation.get_sensitivity_data(nonBaselineScenario, "food_demand_prices",
                                                                          SSP, RCP=RCP, source="masked")
    released_staple_expenditure = released_staple_expenditure[
        released_staple_expenditure[['input']].isin(["FoodDemand_Staples"]).any(axis=1)]
    pyrolysis_staple_expenditure = pyrolysis_staple_expenditure[
        pyrolysis_staple_expenditure[['input']].isin(["FoodDemand_Staples"]).any(axis=1)]

    # get food consumption
    released_staple_consumption = data_manipulation.get_sensitivity_data(["released"], "food_demand_per_capita", SSP,
                                                                         RCP=RCP, source="original")
    pyrolysis_staple_consumption = data_manipulation.get_sensitivity_data(nonBaselineScenario, "food_demand_per_capita",
                                                                          SSP, RCP=RCP, source="masked")
    released_staple_consumption = released_staple_consumption[
        released_staple_consumption[['input']].isin(["FoodDemand_Staples"]).any(axis=1)]
    pyrolysis_staple_consumption = pyrolysis_staple_consumption[
        pyrolysis_staple_consumption[['input']].isin(["FoodDemand_Staples"]).any(axis=1)]

    # get GDP per capita
    released_GDP_capita = data_manipulation.get_sensitivity_data(nonBaselineScenario, "GDP_per_capita_PPP_by_region",
                                                                 SSP, RCP=RCP, source="original")
    pyrolysis_GDP_capita = data_manipulation.get_sensitivity_data(nonBaselineScenario, "GDP_per_capita_PPP_by_region",
                                                                  SSP, RCP=RCP, source="masked")

    # calculate consumption times price divided by GDP per capita
    released_consumption = pd.merge(released_staple_consumption, released_staple_expenditure, how="inner",
                                    on=["SSP", "GCAM", "technology"],
                                    suffixes=("_pcal", "_$"))
    pyrolysis_consumption = pd.merge(pyrolysis_staple_consumption, pyrolysis_staple_expenditure, how="inner",
                                     on=["SSP", "GCAM", "technology"],
                                     suffixes=("_pcal", "_$"))

    # other scaling factors
    released_FA = pd.merge(released_consumption, released_GDP_capita, how="left", on=["SSP", "GCAM"],
                           suffixes=("", "_capita"))
    pyrolysis_FA = pd.merge(pyrolysis_consumption, pyrolysis_GDP_capita, how="left", on=["SSP", "GCAM"],
                            suffixes=("", "_capita"))
    released_FA = pd.merge(released_FA, released_caloric_consumption, how="left", on=["SSP", "GCAM"],
                           suffixes=("", "_caloric"))
    pyrolysis_FA = pd.merge(pyrolysis_FA, pyrolysis_caloric_consumption, how="left", on=["SSP", "GCAM"],
                            suffixes=("", "_caloric"))
    for i in c.GCAMConstants.biochar_x:
        # released_FA[str(i)] corresponds to the capita column
        released_FA[str(i)] = released_consumption[str(i) + "_pcal"] * released_consumption[str(i) + "_$"] / \
                              released_FA[str(i)] * 3.542 / released_FA[str(i) + "_caloric"]
        pyrolysis_FA[str(i)] = pyrolysis_consumption[str(i) + "_pcal"] * pyrolysis_consumption[str(i) + "_$"] / \
                               pyrolysis_FA[str(i)] * 3.542 / pyrolysis_FA[str(i) + "_caloric"]

    perc_diff_FA = data_manipulation.percent_difference(released_FA, pyrolysis_FA, ["GCAM", "SSP"])
    plotting.plot_world(perc_diff_FA, ["%"], SSP, "year", "Units", c.GCAMConstants.biochar_x,
                        "Food Accessibility near midcentury ", RCP, nonBaselineScenario)

    # calculate pcal per capita
    released_pcal_pop = pd.merge(released_Pcal, released_pop, how="inner", on=["SSP", "GCAM"],
                                 suffixes=("_pcal", "_pop"))
    pyrolysis_pcal_pop = pd.merge(pyrolysis_Pcal, pyrolysis_pop, how="inner", on=["SSP", "GCAM"],
                                  suffixes=("_pcal", "_pop"))

    # calculate pcal per capita in the biochar year
    released_pcal_pop["pcal_capita_" + biochar_year] = released_pcal_pop[biochar_year + "_pcal"] / (
                1000 * released_pcal_pop[
            biochar_year + "_pop"]) * 1000000000000 / 365 / 2  # * peta to kilo/365/conversion factor of 2 randomly
    pyrolysis_pcal_pop["pcal_capita_" + biochar_year] = pyrolysis_pcal_pop[biochar_year + "_pcal"] / (
            1000 * pyrolysis_pcal_pop[biochar_year + "_pop"]) * 1000000000000 / 365 / 2
    released_pcal_pop["Units"] = "kcal/capita/day"
    pyrolysis_pcal_pop["Units"] = "kcal/capita/day"

    merged_pcal = released_pcal_pop.merge(pyrolysis_pcal_pop, how="inner", on=["SSP", "GCAM", "technology_pcal"],
                                          suffixes=("_left", "_right"))
    merged_pcal["pcal_capita_" + biochar_year] = merged_pcal["pcal_capita_" + biochar_year + "_right"] - merged_pcal[
        "pcal_capita_" + biochar_year + "_left"]
    data_manipulation.drop_missing(merged_pcal).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(
            RCP) + "/change_in_consumption_kcal_capita_day.csv")

    # extract population and identifying information in biochar_year for weighted average calculations
    merged_pop = pd.DataFrame()
    merged_pop["pcal_capita_" + biochar_year] = merged_pcal[biochar_year + "_pop_right"]
    merged_pop["GCAM"] = merged_pcal["GCAM"]
    merged_pop["SSP"] = merged_pcal["SSP"]
    merged_pop["technology_pcal"] = merged_pcal["technology_pcal"]

    plotting.plot_regional_vertical_avg(merged_pcal, "pcal_capita_" + biochar_year, SSP,
                                        "change in food demand (kcal/person/day)",
                                        "change in food demand in " + biochar_year + " in " + str(SSP[0]),
                                        "technology_pcal", merged_pop, RCP, nonBaselineScenario)




def luc_by_region(nonBaselineScenario, RCP, SSP, biochar_year):
    """
    plots information related to land use changes between pyrolysis and reference scenario
    :param nonBaselineScenario: the scenario to be compared to the released scenario
    :param RCP: the RCP pathways being considered
    :param SSP: the SSP pathways being considered
    :param biochar_year: the year the biochar pathways are being evaluated
    :return: N/A
    """
    # get luc data
    released_luc = data_manipulation.get_sensitivity_data(["released"], "LUC_emissions_by_LUT", SSP, RCP=RCP,
                                                          source="original")
    pyrolysis_luc = data_manipulation.get_sensitivity_data(nonBaselineScenario, "LUC_emissions_by_LUT", SSP, RCP=RCP,
                                                           source="masked")

    released_luc = data_manipulation.group(released_luc, ["GCAM", "SSP"])
    pyrolysis_luc = data_manipulation.group(pyrolysis_luc, ["GCAM", "SSP"])
    flat_diff_luc = data_manipulation.flat_difference(released_luc, pyrolysis_luc, ["GCAM", "SSP"])
    perc_diff_luc = data_manipulation.percent_difference(released_luc, pyrolysis_luc, ["GCAM", "SSP"])
    data_manipulation.drop_missing(flat_diff_luc.drop("LandLeaf", axis=1)).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(RCP) + "/change_in_LUC_emissions.csv")
    data_manipulation.drop_missing(perc_diff_luc.drop("LandLeaf", axis=1)).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(RCP) + "/percent_change_in_LUC_emissions.csv")

    plotting.plot_world_by_years(flat_diff_luc, ["MtC/yr"], "Units", c.GCAMConstants.biochar_x, SSP,
                                 "net difference in LUC emissions by region", RCP, nonBaselineScenario)

    flat_diff_luc = data_manipulation.group(flat_diff_luc, ["SSP"])
    plotting.plot_line_by_product(flat_diff_luc, ["SSP1"], "SSP", ["SSP1"], "SSP",
                                  "Net LUC compared to reference scenario", RCP, nonBaselineScenario)

    released_luc = data_manipulation.get_sensitivity_data(["released"], "LUC_emissions_by_LUT", SSP, RCP=RCP,
                                                          source="original")
    pyrolysis_luc = data_manipulation.get_sensitivity_data(nonBaselineScenario, "LUC_emissions_by_LUT", SSP, RCP=RCP,
                                                           source="masked")
    released_luc["LandLeaf"] = released_luc.apply(lambda row: data_manipulation.relabel_detailed_land_use(row), axis=1)
    pyrolysis_luc["LandLeaf"] = pyrolysis_luc.apply(lambda row: data_manipulation.relabel_detailed_land_use(row),
                                                    axis=1)
    released_luc = data_manipulation.group(released_luc, ["GCAM", "SSP", "LandLeaf"])
    pyrolysis_luc = data_manipulation.group(pyrolysis_luc, ["GCAM", "SSP", "LandLeaf"])
    flat_diff_luc = data_manipulation.flat_difference(released_luc, pyrolysis_luc, ["GCAM", "SSP", "LandLeaf"])
    for i in c.GCAMConstants.biochar_x:
        plotting.plot_regional_hist_avg(flat_diff_luc, str(i), SSP, "count region-LandLeaf",
                                        "Flat diffference in LUC emissions between pyrolysis and reference scenario in " + str(i),
                                        "LandLeaf", "na", RCP, nonBaselineScenario)



def pyrolysis_costing(nonBaselineScenario, RCP, SSP, biochar_year):
    """
    returns information related to the cost of the pyrolysis scenario
    :param nonBaselineScenario: the scenario to be compared to the released scenario
    :param RCP: the RCP pathways being considered
    :param SSP: the SSP pathways being considered
    :param biochar_year: the year being evaluated
    :return: N/A
    """
    # get total costs
    total_cost = data_manipulation.get_sensitivity_data(nonBaselineScenario, "costs_by_tech", SSP, RCP=RCP,
                                                        source="masked")
    unit_cost = data_manipulation.get_sensitivity_data(nonBaselineScenario, "costs_by_tech_and_input", SSP, RCP=RCP,
                                                       source="masked")
    feedstock_cost = data_manipulation.get_sensitivity_data(nonBaselineScenario, "prices_of_all_markets", SSP, RCP=RCP,
                                                            source="masked")

    total_cost = total_cost[total_cost[['sector']].isin(['biochar']).any(axis=1)]
    unit_cost = unit_cost[unit_cost[['sector']].isin(['biochar']).any(axis=1)]
    feedstock_cost = feedstock_cost[feedstock_cost[['product']].isin(
        ['beef manure', 'dairy manure', 'goat manure', 'pork manure', 'poultry manure', "biochar"]).any(axis=1)]
    data_manipulation.drop_missing(total_cost[["GCAM", biochar_year, "technology", "Units"]]).to_csv("data/data_analysis/total_cost_pyrolysis.csv")
    data_manipulation.drop_missing(unit_cost[["GCAM", biochar_year, "technology", "Units"]]).to_csv("data/data_analysis/unit_cost_pyrolysis.csv")
    data_manipulation.drop_missing(feedstock_cost[["GCAM", biochar_year, "product", "Units"]]).to_csv("data/data_analysis/feedstock_cost_pyrolysis.csv")

    feedstock_cost = feedstock_cost[feedstock_cost[['product']].isin(
        ['beef manure', 'dairy manure', 'goat manure', 'pork manure', 'poultry manure']).any(axis=1)]
    # drop outliers
    feedstock_cost = feedstock_cost[feedstock_cost[biochar_year] < 3]
    feedstock_cost[biochar_year] = feedstock_cost[biochar_year] / 0.17 * 1000
    feedstock_cost["Units"] = "USD$/ton"
    plotting.plot_regional_hist_avg(feedstock_cost, biochar_year, SSP, "count", "price distribution of manures in " + biochar_year, "product",
                                    "na", RCP, nonBaselineScenario)


def biochar_rate_by_land_size(nonBaselineScenario, RCP, SSP):
    """
    scatter plot of biochar application rate to the size of the biochar land area in the biochar year
    :param nonBaselineScenario: the scenario to be compared to the released scenario
    :param RCP: the RCP pathways being considered
    :param SSP: the SSP pathways being considered
    :return: N/A
    """
    # read in biochar application rates, and get the application rates
    biochar_app_rate = pd.read_csv("gcam/input/gcamdata/inst/extdata/aglu/A_ag_kgbioha_R_C_Y_GLU_irr_level.csv")
    region_names = pd.read_csv("gcam/input/gcamdata/inst/extdata/water/basin_to_country_mapping.csv", skiprows=7)

    # add extra data to dataframe to help downstream code
    biochar_app_rate['GCAM'] = biochar_app_rate['region']
    biochar_app_rate['Units'] = 'kg biochar/ha/yr'
    biochar_app_rate["SSP"] = SSP[0]

    # rename GLU for mapping
    biochar_app_rate = biochar_app_rate.merge(region_names, "left", left_on="GLU", right_on="GLU_code")
    biochar_app_rate = biochar_app_rate[
        ["kg_bio_ha", "Units", "SSP", "region", "GCAM_commodity", "GCAM_subsector", "GLU_name", "Irr_Rfd"]]
    biochar_app_rate["Irr_Rfd"] = biochar_app_rate["Irr_Rfd"].str.upper()

    # extract information on crops
    biochar_app_rate['technology'] = biochar_app_rate['GCAM_commodity']
    biochar_app_rate['GCAM'] = biochar_app_rate['region']
    biochar_app_rate['technology'] = biochar_app_rate.apply(
        lambda row: data_manipulation.relabel_food(row, "technology"), axis=1)

    # read in detailed land allocation
    # biochar cropland application changes
    land_use = data_manipulation.get_sensitivity_data(nonBaselineScenario, "detailed_land_allocation", SSP, RCP=RCP,
                                                      source="masked")
    # get biochar land use type information
    land_use[["GCAM_subsector", "GLU_name", "Irr_Rfd", "MGMT"]] = land_use['LandLeaf'].str.split("_", expand=True)
    land_use = land_use[
        ["GCAM", "GCAM_subsector", "GLU_name", "Irr_Rfd", "MGMT"] + [str(i) for i in c.GCAMConstants.future_x]]
    land_use["Irr_Rfd"] = land_use["Irr_Rfd"].str.upper()

    # merge datasets
    scatter_data = pd.merge(biochar_app_rate, land_use, "left", on=["GCAM", "GCAM_subsector", "GLU_name", "Irr_Rfd"])

    # remove high outlier
    outlier_cutoff = 3000  # kg/ha/yr
    scatter_data = scatter_data[scatter_data['kg_bio_ha'] < outlier_cutoff]

    # keep only biochar land management and relabel crops
    scatter_data = scatter_data[scatter_data[['MGMT']].isin(["biochar"]).any(axis=1)]
    scatter_data["GCAM_subsector"] = scatter_data.apply(
        lambda row: data_manipulation.relabel_land_crops(row, "GCAM_subsector"), axis=1)

    # plot datasets
    for i in c.GCAMConstants.biochar_x:
        # check to ensure that biochar land doesn't exist until biochar is adopted in 2035
        plotting.plot_regional_vertical(scatter_data, str(i), SSP, y_label="land area (thousand km2)",
                                        title="distribution of usage of biochar lands in " + str(i),
                                        x_column="kg_bio_ha",
                                        x_label="kg biochar/ha/yr", y_column="GCAM_subsector", RCP=RCP,
                                        nonBaselineScenario=nonBaselineScenario)


def farmer_economics(nonBaselineScenario, RCP, SSP, biochar_year):
    """
    plots information related to the changes in farming due to biochar production
    :param nonBaselineScenario: the scenario to be compared to the released scenario
    :param RCP: the RCP pathways being considered
    :param SSP: the SSP pathways being considered
    :param biochar_year: the year the biochar pathways are being evaluated
    :return: N/A
    """
    # get data
    pyrolysis_yields = data_manipulation.get_sensitivity_data(nonBaselineScenario, "ag_tech_yield", SSP, RCP=RCP,
                                                              source="masked")
    pyrolysis_land = data_manipulation.get_sensitivity_data(nonBaselineScenario, "detailed_land_allocation", SSP,
                                                            RCP=RCP, source="masked")
    pyrolysis_profit_rate = data_manipulation.get_sensitivity_data(nonBaselineScenario, "profit_rate", SSP, RCP=RCP,
                                                                   source="masked")
    released_yields = data_manipulation.get_sensitivity_data(["released"], "ag_tech_yield", SSP, RCP=RCP,
                                                             source="original")
    released_land = data_manipulation.get_sensitivity_data(["released"], "detailed_land_allocation", SSP, RCP=RCP,
                                                           source="original")
    released_profit_rate = data_manipulation.get_sensitivity_data(["released"], "profit_rate", SSP, RCP=RCP,
                                                                  source="original")
    pyrolysis_yields["Units"] = "NA"
    released_yields["Units"] = "NA"

    # profit rates
    # change in profit to the farmer compared to baseline (hi mgmt type)
    units = pyrolysis_profit_rate["Units"].unique()[0]
    pyrolysis_profit_rate[["Crop", "basin", "rainfed", "mgmt"]] = pyrolysis_profit_rate['LandLeaf'].str.split('_',
                                                                                                              expand=True)
    released_profit_rate[["Crop", "basin", "rainfed", "mgmt"]] = released_profit_rate['LandLeaf'].str.split('_',
                                                                                                            expand=True)
    released_hi_profit = released_profit_rate[released_profit_rate[['mgmt']].isin(["hi"]).any(axis=1)].copy(
        deep=True)
    pyrolysis_biochar_profit = pyrolysis_profit_rate[
        pyrolysis_profit_rate[['mgmt']].isin(["biochar"]).any(axis=1)].copy(
        deep=True)
    pyrolysis_diff_profit = pd.merge(released_hi_profit, pyrolysis_biochar_profit, how="left",
                                     on=["GCAM", "SSP", "Crop", "basin", "rainfed"], suffixes=("_left", "_right"))
    pyrolysis_diff_profit["Units"] = units
    pyrolysis_diff_profit['Crop'] = pyrolysis_diff_profit.apply(
        lambda row: data_manipulation.relabel_land_crops(row, "Crop"), axis=1)
    for i in c.GCAMConstants.x:
        pyrolysis_diff_profit[str(i)] = pyrolysis_diff_profit[str(i) + "_right"] - pyrolysis_diff_profit[
            str(i) + "_left"]

    pyrolysis_diff_profit = pyrolysis_diff_profit.sort_values(by=biochar_year)
    flat_diff_small = pyrolysis_diff_profit[
        (-6e7 < pyrolysis_diff_profit[biochar_year]) & (pyrolysis_diff_profit[biochar_year] < 6e7)]
    flat_diff_large = pyrolysis_diff_profit[
        (-6e7 >= pyrolysis_diff_profit[biochar_year]) | (pyrolysis_diff_profit[biochar_year] >= 6e7)]
    plotting.plot_regional_hist_avg(flat_diff_small, biochar_year, ["SSP1"], "count crop-basin-irrigation-year",
                                    "histogram of small farmer profit rate changes at the crop level in " + biochar_year, "Crop", "na",
                                    RCP, nonBaselineScenario)
    plotting.plot_regional_hist_avg(flat_diff_large, biochar_year, ["SSP1"], "count crop-basin-irrigation-year",
                                    "histogram of large farmer profit rate changes at the crop level in " + biochar_year, "Crop", "na",
                                    RCP, nonBaselineScenario)

    # change in per crop supply
    pyrolysis_yields_lands = pd.merge(pyrolysis_yields, pyrolysis_land, "left", left_on=["GCAM", "SSP", "technology"],
                                      right_on=["GCAM", "SSP", "LandLeaf"], suffixes=("_left", "_right"))
    released_yields_lands = pd.merge(released_yields, released_land, "left",
                                     left_on=["GCAM", "SSP", "technology"], right_on=["GCAM", "SSP", "LandLeaf"],
                                     suffixes=("_left", "_right"))
    pyrolysis_lands_grouping = pyrolysis_yields_lands.copy(deep=True)
    released_lands_grouping = released_yields_lands.copy(deep=True)

    for i in c.GCAMConstants.x:
        pyrolysis_yields_lands[str(i)] = pyrolysis_yields_lands[str(i) + "_left"] * pyrolysis_yields_lands[
            str(i) + "_right"]
        released_yields_lands[str(i)] = released_yields_lands[str(i) + "_left"] * released_yields_lands[
            str(i) + "_right"]
        pyrolysis_lands_grouping[str(i)] = pyrolysis_yields_lands[str(i) + "_right"]
        released_lands_grouping[str(i)] = released_yields_lands[str(i) + "_right"]

    # group by crop
    pyrolysis_yields_lands[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"
    released_yields_lands[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"
    pyrolysis_lands_grouping[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"
    released_lands_grouping[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"
    pyrolysis_yields_lands["sector"] = pyrolysis_yields_lands["sector_left"]
    released_yields_lands["sector"] = released_yields_lands["sector_left"]
    pyrolysis_lands_grouping["sector"] = pyrolysis_lands_grouping["sector_left"]
    released_lands_grouping["sector"] = released_lands_grouping["sector_left"]
    pyrolysis_effective_yield = data_manipulation.group(pyrolysis_yields_lands, ["GCAM", "SSP", "sector"])
    released_effective_yield = data_manipulation.group(released_yields_lands, ["GCAM", "SSP", "sector"])
    pyrolysis_lands_grouping = data_manipulation.group(pyrolysis_lands_grouping, ["GCAM", "SSP", "sector"])
    released_lands_grouping = data_manipulation.group(released_lands_grouping, ["GCAM", "SSP", "sector"])

    # divide by available crop land
    pyrolysis_effective_yield = pd.merge(pyrolysis_effective_yield, pyrolysis_lands_grouping, "left",
                                         on=["GCAM", "SSP", "sector"],
                                         suffixes=("_left", "_right"))
    released_effective_yield = pd.merge(released_effective_yield, released_lands_grouping, "left",
                                        on=["GCAM", "SSP", "sector"],
                                        suffixes=("_left", "_right"))
    for i in c.GCAMConstants.x:
        pyrolysis_effective_yield[str(i)] = pyrolysis_effective_yield[str(i) + "_left"] / pyrolysis_effective_yield[
            str(i) + "_right"]
        released_effective_yield[str(i)] = released_effective_yield[str(i) + "_left"] / released_effective_yield[
            str(i) + "_right"]

    pyrolysis_effective_yield[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"
    released_effective_yield[
        ['Version', 'output', 'concentration', 'input', 'product', 'fuel', 'LandLeaf', 'GHG', "Units", "subsector",
         "technology"]] = "NA"

    pyrolysis_effective_yield = pyrolysis_effective_yield[c.GCAMConstants.column_order]
    released_effective_yield = released_effective_yield[c.GCAMConstants.column_order]

    # yield differences between crops
    flat_diff_effective_yields = data_manipulation.flat_difference(released_effective_yield, pyrolysis_effective_yield,
                                                                   ["GCAM", "SSP", "sector"])
    plotting.plot_regional_hist_avg(flat_diff_effective_yields, biochar_year, ["SSP1"], "count region-year",
                                    "histogram of yield changes at the crop level in " + biochar_year, "sector", "na", RCP,
                                    nonBaselineScenario)

    empty_pyro = pyrolysis_land[pyrolysis_land["LandLeaf"].str.contains("biochar")].copy(deep=True)
    for i in c.GCAMConstants.x:
        empty_pyro[str(i)] = 0
    released_land = pd.concat([released_land, empty_pyro])
    # subtract from pyrolysis land
    flat_diff_mgmt = data_manipulation.flat_difference(released_land, pyrolysis_land, ["GCAM", "SSP", "LandLeaf"])
    flat_diff_mgmt[["Crop", "basin", "rainfed", "mgmt"]] = flat_diff_mgmt['LandLeaf'].str.split('_', expand=True)
    flat_diff_mgmt = flat_diff_mgmt[flat_diff_mgmt['mgmt'].notna()]
    flat_diff_mgmt = flat_diff_mgmt.sort_values(by=biochar_year)
    flat_diff_small = flat_diff_mgmt[(-1 < flat_diff_mgmt[biochar_year]) & (flat_diff_mgmt[biochar_year] < 1)]
    flat_diff_large = flat_diff_mgmt[(-1 >= flat_diff_mgmt[biochar_year]) | (flat_diff_mgmt[biochar_year] >= 1)]
    plotting.plot_regional_hist_avg(flat_diff_small, biochar_year, ["SSP1"], "count basin-crop-irrigation",
                                    "histogram of small land mgmt changes in terms of area compared to reference scenario in " + biochar_year,
                                    "mgmt", "na", RCP, nonBaselineScenario)
    plotting.plot_regional_hist_avg(flat_diff_large, biochar_year, ["SSP1"], "count basin-crop-irrigation",
                                    "histogram of large land mgmt changes in terms of area compared to reference scenario in " + biochar_year,
                                    "mgmt", "na", RCP, nonBaselineScenario)

    # land leaf shares histogram
    pyrolysis_landleafs = data_manipulation.get_sensitivity_data(nonBaselineScenario, "land_leaf_shares", SSP, RCP=RCP,
                                                                 source="masked")
    pyrolysis_landleafs[["Crop", "basin", "rainfed", "MGMT"]] = pyrolysis_landleafs['LandLeaf'].str.split('_',
                                                                                                          expand=True)
    pyrolysis_landleafs = pyrolysis_landleafs[pyrolysis_landleafs[['MGMT']].isin(["biochar"]).any(axis=1)]
    pyrolysis_landleafs["Units"] = "land share"
    for i in c.GCAMConstants.biochar_x:
        plotting.plot_regional_hist_avg(pyrolysis_landleafs, str(i), SSP, "count land leafs",
                                        "histogram of land leaf shares for biochar lands in " + str(i), "Crop", "na",
                                        RCP, nonBaselineScenario)

    # output differences in carbon prices
    c_pyro_price = data_manipulation.get_sensitivity_data(nonBaselineScenario, "CO2_prices", SSP, RCP=RCP,
                                                          source="masked")
    c_rel_price = data_manipulation.get_sensitivity_data(["released"], "CO2_prices", SSP, RCP=RCP,
                                                         source="original")
    product = ["CO2"]
    c_rel_price = c_rel_price[c_rel_price[['product']].isin(product).any(axis=1)]
    c_pyro_price = c_pyro_price[c_pyro_price[['product']].isin(product).any(axis=1)]
    for i in c.GCAMConstants.x:
        c_rel_price[str(i)] = c_rel_price[
                                  str(i)] * 2.42  # https://data.bls.gov/cgi-bin/cpicalc.pl?cost1=1.00&year1=199001&year2=202401
        c_pyro_price[str(i)] = c_pyro_price[
                                   str(i)] * 2.42  # https://data.bls.gov/cgi-bin/cpicalc.pl?cost1=1.00&year1=199001&year2=202401
    flat_diff_c_price = data_manipulation.flat_difference(c_pyro_price, c_rel_price, ["SSP", "GCAM"])
    flat_diff_c_price["Units"] = "USD$2024/t C"
    perc_diff_c_price = data_manipulation.percent_difference(c_pyro_price, c_rel_price, ["SSP", "GCAM"])
    data_manipulation.drop_missing(flat_diff_c_price).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(
            RCP) + "/change_in_carbon_price.csv")
    data_manipulation.drop_missing(perc_diff_c_price).to_csv(
        "data/data_analysis/supplementary_tables/" + str(nonBaselineScenario) + "/" + str(
            RCP) + "/percent_change_in_carbon_price.csv")


def main():
    """
    Main method for running all scripts
    :return: N/A
    """
    reference_SSP = ["SSP1"]
    reference_RCP = "baseline"
    other_scenario = ["Baseline",
                      "HighBiocharCost",
                      "HighBiocharNUE",
                      "HighBiocharNutrients",
                      "HighBiocharSoilN2O",
                      "HighBiocharYield",
                      "HighCropYield",
                      "HighGCAMLandShare",
                      "HighGCAMManurePrice",
                      "LowBiocharNutrients",
                      "LowBiocharCost",
                      "LowBiocharNUE",
                      "LowBiocharSoilN2O",
                      "LowBiocharYield",
                      "LowCropYield",
                      "LowGCAMLandShare",
                      "LowGCAMManurePrice",
                      "Highadoption70",
                      "HighCarbonStability",
                      "Lowadoption30",
                      "LowCarbonStability"]
    biochar_year = "2050"
    #biochar_rate_by_land_size(other_scenario, reference_RCP, reference_SSP)
    farmer_economics(other_scenario, reference_RCP, reference_SSP, biochar_year)
    #pyrolysis_costing(other_scenario, reference_RCP, reference_SSP, biochar_year)
    #animal_feed_and_products(other_scenario, reference_RCP, reference_SSP, biochar_year)
    #luc_by_region(other_scenario, reference_RCP, reference_SSP, biochar_year)
    pop_and_calories(other_scenario, reference_RCP, reference_SSP, biochar_year)


if __name__ == '__main__':
    main()
