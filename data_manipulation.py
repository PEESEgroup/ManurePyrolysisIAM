import numpy as np
import constants as c
import pandas as pd
import scipy.stats as stats


def flat_difference(old, new, columns):
    """
    Calculates the flat difference between two dataframes (new - old)
    :param old: the old dataframe
    :param new: the new dataframe
    :param columns: the list of columns that will uniquely identify each product
    :return: a combined dataframe
    """
    # get values from dataframes
    merged = old.merge(new, how="left", on=columns, suffixes=("_left", "_right"))

    for i in c.GCAMConstants.future_x:
        merged[str(i)] = merged[str(i) + "_right"] - merged[str(i) + "_left"]
        merged = merged.drop([str(i) + "_left"], axis=1)
        merged = merged.drop([str(i) + "_right"], axis=1)

    # update the version name
    # merged['Comparison'] = "change between " + merged.apply(lambda row: str_version(row), axis=1)

    # fix the column names
    merged.columns = merged.columns.str.replace("_right", '')  # only label columns
    merged = merged[c.GCAMConstants.column_order]

    return merged


def percent_difference(old, new, columns):
    """
    calculates the percent difference between two dataframes
    :param old: the old dataframe
    :param new: the new dataframe
    :param columns: the columns necessary to uniquely identify a product
    :return: a dataframe containing the percent change between the dataframes
    """
    # get values from dataframes
    merged = old.merge(new, how="left", on=columns, suffixes=("_left", "_right"))

    for i in c.GCAMConstants.future_x:
        merged[str(i)] = merged.apply(lambda row: calc_percs(row, i), axis=1)
        merged = merged.drop([str(i) + "_left"], axis=1)
        merged = merged.drop([str(i) + "_right"], axis=1)

    # update columns
    merged = merged.drop(['Units_left'], axis=1)
    merged = merged.drop(['Units_right'], axis=1)
    merged['Units'] = '%'
    # merged['Comparison'] = "% change between " + merged.apply(lambda row: str_version(row), axis=1)

    # replace columns
    merged.columns = merged.columns.str.replace("_right", '')  # only label columns
    merged = merged[c.GCAMConstants.column_order]

    return merged


def str_version(row):
    return str(row['Version_right']) + " and " + str(row['Version_left'])


def calc_percs(row, i):
    """
    if a row has no baseline value, return np.nan instead of percentage change
    :param row: row in dataframe
    :param i: string of product year
    :return: percent difference from right entry over left
    """
    if row[str(i) + "_left"] == 0:
        return np.nan
    else:
        return 100 * (row[str(i) + "_right"] - row[str(i) + "_left"]) / (row[str(i) + "_left"] + 1e-7)


def percent_of_total(old, new, columns, biochar_year):
    """
    calculates the percent difference between two dataframes
    :param old: the old dataframe
    :param new: the new dataframe
    :param columns: the columns necessary to uniquely identify a product
    :param biochar_year: a year with widespread biochar adoption
    :return: a dataframe containing the percent change between the dataframes
    """
    # get values from dataframes
    merged = old.merge(new, how="left", on=columns, suffixes=("_left", "_right"))
    df = pd.DataFrame()

    for j in merged["Version_right"].unique().tolist():
        for k in merged["GCAM"].unique().tolist():
            merged_filter = merged[(merged['Version_right'].isin([j])) & (merged['GCAM'].isin([k]))]  # filter by region
            for i in c.GCAMConstants.future_x:
                sum_col = merged_filter[str(i) + "_left"].sum()
                merged_filter = merged_filter.assign(
                    e=100 * (merged_filter[str(i) + "_right"] - merged_filter[str(i) + "_left"]) / sum_col)
                merged_filter[str(i)] = merged_filter["e"]
                merged_filter = merged_filter.drop([str(i) + "_left"], axis=1)
                merged_filter = merged_filter.drop([str(i) + "_right"], axis=1)

                if k == "Global" and str(i) == biochar_year:
                    print("total is", sum_col, "thousand sq km in", j)
                    merged_filter["GCAM"] = "Global (net)"

            merged_filter.drop(["e"], axis=1)
            df = pd.concat([df, merged_filter])

    # update columns
    df = df.drop(['Units_left'], axis=1)
    df = df.drop(['Units_right'], axis=1)
    df['Units'] = '%'

    # replace columns
    df.columns = df.columns.str.replace("_right", '')
    df = df[c.GCAMConstants.column_order]

    return df


def group(df, columns):
    """
    Groups a dataframe with many subproducts into a single line via summation
    :param df: the dataframe being grouped
    :param columns: the list of columns used to form a group
    :return: a dataframe with grouped entries
    """
    unit = df.groupby(columns).first().reset_index()["Units"]
    df = df.groupby(columns).sum(min_count=1)
    df = df.reset_index()
    df['Units'] = unit
    df['SSP'] = df.apply(lambda row: relabel_SSP(row), axis=1)
    df = df[c.GCAMConstants.column_order]
    return df


def relabel_SSP(row):
    """
    lambda function to relabel GCAM SSP after grouping.
    :param row: row of data
    :return: ungrouped GCAM SSP
    """
    SSP = row["SSP"]  # Products in staples vs. non staples in A_demand_subsector.csv
    if "1" in SSP:
        return "SSP1"
    elif "2" in SSP:
        return "SSP2"
    elif "3" in SSP:
        return "SSP3"
    elif "4" in SSP:
        return "SSP4"
    elif "5" in SSP:
        return "SSP5"
    else:
        return "error"



def label_sequestration_sectors(row):
    """
    returns the aggregated sector to better show plots
    :param row: row of a pandas dataframe
    :return: aggregated sector or original one
    """
    if row['sector'] in ["H2 central production", "H2 wholesale dispensing"]:
        return "hydrogen"
    elif row['sector'] in ["backup_electricity", "district heat", "refining"]:
        return "other energy sector"
    elif row['sector'] in ["chemical energy use", "chemical feedstocks", "N fertilizer", "alumina", "cement",
                           "iron and steel", "other industrial energy use", "other industrial feedstocks",
                           "process heat cement"]:
        return "industrial energy use"
    elif row['sector'] in ["comm cooling", "comm heating", "comm others"]:
        return "commercial energy use"
    elif row['sector'] in ["elec_biomass (IGCC CCS)", "elec_biomass (IGCC)", "elec_biomass (conv CCS)",
                           "elec_biomass (conv)"]:
        return "electricity - biomass"
    elif row['sector'] in ["elec_coal (IGCC CCS)", "elec_coal (IGCC)", "elec_coal (conv pul CCS)",
                           "elec_coal (conv pul)"]:
        return "electricity - coal"
    elif row['sector'] in ["elec_gas (CC CCS)", "elec_gas (CC)", "elec_gas (steam/CT)"]:
        return "electricity - gas"
    elif row['sector'] in ["elec_refined liquids (CC CCS)", "elec_refined liquids (CC)",
                           "elec_refined liquids (steam/CT)"]:
        return "electricity - refined liquids"
    elif row['sector'] in ["gas pipeline", "gas processing"]:
        return "gas processing"
    elif row['sector'] in ["resid cooling", "resid heating", "resid others"]:
        return "commercial energy use"
    elif row['sector'] in ["trn_aviation_intl", "trn_freigh", "trn_freight_road", "trn_pass", "trn_pass_road",
                           "trn_pass_road_LDV", "trn_pass_road_LDV_4W", "trn_shipping_intl", "trn_freight"]:
        return "transportation"
    elif row['sector'] in ['regional biomass', 'regional biomassOil', "regional corn for ethanol",
                           "regional sugar for ethanol"]:
        return "other biomass for refining"
    else:
        return row['sector']


def remove__(row, column):
    """
    relabels similar technologies to enable grouping
    :param row: a pd Series from a dataframe
    :param column: the column of the pd series being searched
    :return: the relabeled technology
    """
    to_return = row[column].replace("_", " ")
    return to_return


def relabel_region(row):
    """
    lambda function to relabel GCAM regions for greater accessibility
    :param row: row of data
    :return: updated name of GCAM region
    """
    GCAM_region = row["GCAM"]
    matches = ["Argentina", "Brazil", "Canada", "Central America and Caribbean", "Central Asia", "China",
               "Colombia", "European Free Trade Association", "India", "Indonesia", "Mexico",
               "Japan", "Middle East", "Pakistan", "Russia", "South Africa", "South Asia",
               "Southeast Asia", "South Korea", "Taiwan", "USA", "Global"]
    if any(x in GCAM_region for x in matches):
        return GCAM_region
    elif GCAM_region == "Middle East":
        return "Middle East"
    elif GCAM_region == "Indonesia":
        return "Indonesia"
    elif GCAM_region == "Global":
        return "Global"
    elif "Africa_Eastern" == GCAM_region:
        return "Eastern Africa"
    elif "Africa_Northern" == GCAM_region:
        return "Northern Africa"
    elif "Africa_Southern" == GCAM_region:
        return "Southern Africa"
    elif "Africa_Western" == GCAM_region:
        return "Western Africa"
    elif "Australia_NZ" == GCAM_region:
        return "Australia and New Zealand"
    elif "EU-12" == GCAM_region:
        return "Northeastern EU"
    elif "EU-15" == GCAM_region:
        return "Western EU"
    elif "Europe_Eastern" == GCAM_region:
        return "Eastern Europe"
    elif "Europe_Non_EU" == GCAM_region:
        return "Other Europe"
    elif "South America_Northern" == GCAM_region:
        return "Northern South America"
    elif "South America_Southern" == GCAM_region:
        return "Southern South America"
    elif "Global (net)" == GCAM_region:
        return "Global (net)"
    else:
        return "error"


def relabel_land_use(row):
    """
    lambda function to relabel GCAM LandLeaf for greater accessibility. From: https://jgcri.github.io/gcam-doc/land.html
    :param row: row of data
    :return: updated name of GCAM region
    """
    luc = row["LandLeaf"]
    if luc == "crops":
        return "crops"
    elif luc == "biomass":
        return "biomass for energy"
    elif luc == "grass":
        return "grass land"
    elif luc == "shrubs":
        return "shrub land"
    elif "Hardwood_Forest" == luc or "Softwood_Forest" == luc:
        return "commercial forest"
    elif "UnmanagedHardwood_Forest" == luc or "UnmanagedSoftwood_Forest" == luc:
        return "forest"
    elif "pasture (grazed)" == luc:
        return "intensively-grazed pasture"
    elif "pasture (other)" == luc:
        return "other pasture"
    elif "otherarable" == luc:
        return "other arable land"
    else:
        return luc


def relabel_detailed_land_use(row):
    """
    lambda function to relabel GCAM LandLeaf for greater accessibility. From: https://jgcri.github.io/gcam-doc/land.html
    :param row: row of data
    :return: updated name of GCAM LandLeaf
    """
    luc = row["LandLeaf"]
    if "Grassland" in luc:
        return "Grasslands"
    elif "ProtectedUnmanagedPasture" in luc:
        return "Protected Lands"
    elif "Vegetables" in luc:
        return "Crops"
    elif "FodderHerb" in luc:
        return "Animal Feed"
    elif "MiscCrop" in luc:
        return "Crops"
    elif "OtherGrainC4" in luc:
        return "Crops"
    elif "PalmFruit" in luc:
        return "Crops"
    elif "FiberCrop" in luc:
        return "Crops"
    elif "NutsSeeds" in luc:
        return "Crops"
    elif "OtherGrain" in luc:
        return "Crops"
    elif "Soybean" in luc:
        return "Crops"
    elif "FodderGrass" in luc:
        return "Animal Feed"
    elif "ProtectedGrassland" in luc:
        return "Protected Lands"
    elif "Fruits" in luc:
        return "Crops"
    elif "FodderHerbC4" in luc:
        return "FodderHerb"
    elif "ProtectedUnmanagedForest" in luc:
        return "Protected Lands"
    elif "biomassTree" in luc:
        return "Biomass for Energy"
    elif "OilPalm" in luc:
        return "Crops"
    elif "OtherArableLand" in luc:
        return "Other Arable Land"
    elif "MiscCropTree" in luc:
        return "Crops"
    elif "OilPalmTree" in luc:
        return "Crops"
    elif "Rice" in luc:
        return "Crops"
    elif "Legumes" in luc:
        return "Crops"
    elif "NutsSeedsTree" in luc:
        return "Crops"
    elif "OilCropTree" in luc:
        return "Crops"
    elif "UrbanLand" in luc:
        return "Urban"
    elif "RockIceDesert" in luc:
        return "Rock and Desert"
    elif "RootTuber" in luc:
        return "Crops"
    elif "Corn" in luc:
        return "Crops"
    elif "FruitsTree" in luc:
        return "Crops"
    elif "OilCrop" in luc:
        return "Crops"
    elif "ProtectedShrubland" in luc:
        return "Protected Lands"
    elif "SugarCrop" in luc:
        return "Crops"
    elif "UnmanagedForest" in luc:
        return "Forest"
    elif "SugarCropC4" in luc:
        return "Crops"
    elif "Pasture" in luc:
        return "Pasture"
    elif "Forest" in luc:
        return "Forest"
    elif "biomassGrass" in luc:
        return "Biomass for Energy"
    elif "Shrubland" in luc:
        return "Shrubland"
    elif "UnmanagedPasture" in luc:
        return "Pasture"
    elif "Tundra" in luc:
        return "Tundra"
    elif "Wheat" in luc:
        return "Crops"
    elif "CornC4" in luc:
        return "Crops"
    return "error"


def relabel_food(row, column):
    """
    lambda function to relabel GCAM food categories for greater accessibility. From: A_demand_technology.csv
    :param row: row of data
    :param column: column of data to be processed
    :return: updated name of GCAM region
    """
    food = row[column]
    if food == "NutsSeeds":
        return "Nuts and Seeds"
    elif "Oil" == food:
        return "Plant Oils"
    elif "OilCrop" == food:
        return "Plant Oils"
    elif food == "OilPalm":
        return "Palm Oil"
    elif "FiberCrop" == food:
        return "Fiber Crops"
    elif "MiscCrop" == food:
        return "Other Crops"
    elif "OtherGrain" == food:
        return "Other Grains"
    elif "RootTuber" == food:
        return "Roots and Tubers"
    elif "SheepGoat" == food:
        return "Sheep and Goat"
    elif "OtherMeat_Fish" == food:
        return "Fish"
    elif "FodderGrass" == food:
        return "Fodder Grass"
    elif "FodderHerb" == food:
        return "Fodder Herb"
    elif "SugarCrop" == food:
        return "Sugar Crops"
    elif "FruitsVeg" == food:
        return "Fruits and Vegetables"
    else:
        return food


def relabel_food_demand(row):
    """
    lambda function to relabel GCAM food demand categories for greater accessibility.
    :param row: row of data
    :return: updated name of GCAM region
    """
    food = row["input"]  # Products in staples vs. non staples in A_demand_subsector.csv
    if food == "FoodDemand_NonStaples":
        return "Non-Staples"
    elif food == "FoodDemand_Staples":
        return "Staples"
    else:
        return "error"


def relabel_land_crops(row, column):
    """
    lambda function to relabel GCAM LandLeaf to extract different crop classes
    :param row: row of data
    :param column: column of data
    :return: updated name of GCAM LandLeaf
    """
    luc = row[column]
    if "Grassland" in luc:
        return "grass"
    elif "ProtectedUnmanagedPasture" in luc:
        return "pasture (other)"
    elif "Vegetables" in luc:
        return "Vegetables"
    elif "FodderHerb" in luc:
        return "Animal Feed"
    elif "MiscCrop" in luc:
        return "Other Crops"
    elif "OtherGrainC4" in luc:
        return "Other Crops"
    elif "PalmFruit" in luc:
        return "Fruits"
    elif "FiberCrop" in luc:
        return "Fiber Crops"
    elif "NutsSeeds" in luc:
        return "Nuts and Seeds"
    elif "OtherGrain" in luc:
        return "Other Crops"
    elif "Soybean" in luc:
        return "Soybean"
    elif "FodderGrass" in luc:
        return "Animal Feed"
    elif "ProtectedGrassland" in luc:
        return "grass"
    elif "Fruits" in luc:
        return "Fruits"
    elif "FodderHerbC4" in luc:
        return "Animal Feed"
    elif "ProtectedUnmanagedForest" in luc:
        return "forest (unmanaged)"
    elif "biomassTree" in luc:
        return "Biomass for Energy"
    elif "OilPalm" in luc:
        return "Plant Oils"
    elif "OtherArableLand" in luc:
        return "otherarable"
    elif "MiscCropTree" in luc:
        return "Other Crops"
    elif "OilPalmTree" in luc:
        return "Plant Oils"
    elif "Rice" in luc:
        return "Rice"
    elif "Legumes" in luc:
        return "Legumes"
    elif "NutsSeedsTree" in luc:
        return "Nuts and Seeds"
    elif "OilCropTree" in luc:
        return "Plant Oils"
    elif "UrbanLand" in luc:
        return "urban"
    elif "RockIceDesert" in luc:
        return "rock and desert"
    elif "RootTuber" in luc:
        return "Roots and Tubers"
    elif "Corn" in luc:
        return "Corn"
    elif "FruitsTree" in luc:
        return "Fruits"
    elif "OilCrop" in luc:
        return "Plant Oils"
    elif "ProtectedShrubland" in luc:
        return "shrubs"
    elif "SugarCrop" in luc:
        return "Sugar Crops"
    elif "UnmanagedForest" in luc:
        return "forest (unmanaged)"
    elif "SugarCropC4" in luc:
        return "Sugar Crops"
    elif "Pasture" in luc:
        return "pasture (grazed)"
    elif "Forest" in luc:
        return "forest (managed)"
    elif "biomassGrass" in luc:
        return "Biomass for Energy"
    elif "Shrubland" in luc:
        return "shrubs"
    elif "UnmanagedPasture" in luc:
        return "pasture (other)"
    elif "Tundra" in luc:
        return "tundra"
    elif "Wheat" in luc:
        return "Wheat"
    elif "CornC4" in luc:
        return "Corn"
    return "error"


def relabel_MGMT(row):
    """
    lambda function to relabel GCAM food demand categories for greater accessibility.
    :param row: row of data
    :return: updated name of GCAM region
    """
    item = row["MGMT"]
    if item == "hi":
        return "High Intensity"
    elif item == "lo":
        return "Low Intensity"
    else:
        return "Biochar Application"


def mgmt_to_color(row):
    """
    give colors to the different management types
    :param row: row of a pandas dataframe
    :return: hex code for a color
    """
    if "High Intensity" in row["Management"]:
        return "#698FC6"
    elif "Low Intensity" in row["Management"]:
        return "#BFBE43"
    elif "Biochar Application" in row["Management"]:
        return "#C16861"
    elif "Other Land Use Types" in row["Management"]:
        return "#74A751"


def relabel_region_alluvial(row, biochar_year, base_year):
    """
    process data for the alluvial figure
    :param row: row from a pandas dataframe
    :param biochar_year: year with widespread biochar adoption for analysis
    :param base_year: year with no biochar adoption for analysis
    :return: get data based on the year
    """
    if row["Region_" + str(biochar_year)] is None:
        return row["Region_" + str(base_year)]
    else:
        return row["Region_" + str(biochar_year)]


def relabel_management_alluvial(row, counts, column):
    """
    relabel land management type
    :param row: row of a pandas dataframe
    :param counts: count of the number of entries of a value in the pandas dataframe
    :param column: column in which data is stored
    :return: a string describing the land type and amount
    """
    if row[column] is None:
        return "Other Land Use Types"
    else:
        return row[column] + ":<br>" + f'{float(f"{counts[row.at[column]]:.5g}"):g}' + " km<sup>2</sup>"


def process_luc(land_use, scale_factor, base_year, biochar_year):
    """
    process land use change, attempting to map other land use types to other land use types within a region/basin,
    while matching the same type of crops in both years
    :param land_use: dataframe containing land use types
    :param scale_factor: scale factor increasing the number of rows, and thus increasing the precision of the data,
    while increasing the computation time for the method
    :param base_year: base year for the analysis
    :param biochar_year: year with widespread adoption of biochar for the analysis
    :return: pandas dataframe containing land use mappings for different crops
    """
    land_for_alluvial = pd.DataFrame()
    for r in land_use['GCAM'].unique():
        one_region = land_use[land_use[['GCAM']].isin([r]).any(axis=1)]
        for basin in one_region["Basin"].unique():
            land_per_basin = pd.DataFrame()
            b = one_region[one_region[['Basin']].isin([basin]).any(axis=1)]
            for crop_type in b['Crop'].unique():
                crops = b[b[['Crop']].isin([crop_type]).any(axis=1)]
                by_mgmt_type_2020 = pd.Series()
                by_mgmt_type_2050 = pd.Series()
                for mgmt_type in crops['MGMT'].unique():
                    mgmt = crops[crops[['MGMT']].isin([mgmt_type]).any(axis=1)]
                    year_2020 = pd.Series()
                    year_2050 = pd.Series()
                    regions = mgmt[mgmt[['GCAM']].isin([r]).any(axis=1)]
                    repeat_times_2020 = regions[base_year]
                    repeat_string_2020 = str(crop_type) + "_" + str(mgmt_type) + "_" + str(r)
                    # repeat the land area a number of times equal to the amount of land
                    crop_by_mgmt = pd.Series([repeat_string_2020]).repeat(repeat_times_2020 * scale_factor)
                    # add the repeated number of rows to a list
                    year_2020 = pd.concat([year_2020, crop_by_mgmt])
                    year_2020 = year_2020.reset_index(drop=True)

                    repeat_times_2050 = regions[biochar_year]
                    repeat_string_2050 = str(crop_type) + "_" + str(mgmt_type) + "_" + str(r)
                    # repeat the land area a number of times equal to the amount of land
                    crop_by_mgmt = pd.Series([repeat_string_2050]).repeat(repeat_times_2050 * scale_factor)
                    # add the repeated number of rows to a list
                    year_2050 = pd.concat([year_2050, crop_by_mgmt])
                    year_2050 = year_2050.reset_index(drop=True)

                    # add the lists to a dataframe
                    by_mgmt_type_2020 = pd.concat([by_mgmt_type_2020, year_2020])
                    by_mgmt_type_2050 = pd.concat([by_mgmt_type_2050, year_2050])
                    by_mgmt_type_2050 = by_mgmt_type_2050.reset_index(drop=True)
                    by_mgmt_type_2020 = by_mgmt_type_2020.reset_index(drop=True)

                # gather excess land not present in either year
                by_mgmt_type = pd.DataFrame()
                if by_mgmt_type_2050.size > by_mgmt_type_2020.size:
                    by_mgmt_type[base_year] = by_mgmt_type_2020
                    by_mgmt_type[biochar_year] = by_mgmt_type_2050
                else:
                    by_mgmt_type[base_year] = by_mgmt_type_2020
                    by_mgmt_type[biochar_year] = by_mgmt_type_2050
                by_mgmt_type = by_mgmt_type.fillna("Other Land Use Types")
                by_mgmt_type = by_mgmt_type.reset_index(drop=True)

                land_per_basin = pd.concat([land_per_basin, by_mgmt_type])

            # on a per-basin basis, rearrange luc so that unmanaged land is matched with unmanaged land
            df_both_managed = land_per_basin[
                (land_per_basin[base_year] != "Other Land Use Types") & (land_per_basin[biochar_year] != "Other Land Use Types")].reset_index(drop=True)
            df_2020_managed = land_per_basin[(land_per_basin[biochar_year] == "Other Land Use Types")].reset_index(drop=True)
            df_2050_managed = land_per_basin[(land_per_basin[base_year] == "Other Land Use Types")].reset_index(drop=True)

            # shuffle both lists
            df_2020_managed = df_2020_managed.sample(frac=1, random_state=1).reset_index(drop=True)
            df_2050_managed = df_2050_managed.sample(frac=1, random_state=1).reset_index(drop=True)

            # get excess other land use types and put them in their own df, so that 2020/2050 managed have the same length
            excess = pd.DataFrame()
            if len(df_2020_managed.index) > len(df_2050_managed.index):
                excess = pd.concat([excess, df_2020_managed.iloc[len(df_2050_managed.index):]])
                df_2020_managed = df_2020_managed.iloc[:len(df_2050_managed.index)]
            else:
                excess = pd.concat([excess, df_2050_managed.iloc[len(df_2020_managed.index):]])
                df_2050_managed = df_2050_managed.iloc[:len(df_2020_managed.index)]

            # turn dfs into lists
            df_2020_managed_list = df_2020_managed.to_numpy().tolist()
            df_2050_managed_list = df_2050_managed.to_numpy().tolist()

            pairs = zip(df_2020_managed_list, df_2050_managed_list)
            paired = pd.DataFrame(columns=[base_year, biochar_year])
            for i, j in pairs:
                paired.loc[len(paired)] = [i[0], j[1]]  # based on column order in dataframes

            df_doubly_unmanaged = paired[
                (paired[base_year] == "Other Land Use Types") & (
                            paired[biochar_year] == "Other Land Use Types")].reset_index(drop=True)

            if len(df_doubly_unmanaged) > 0:
                print(r, basin)
            df_both_managed = pd.concat([df_both_managed, paired, excess])
            land_for_alluvial = pd.concat([land_for_alluvial, df_both_managed])

    return land_for_alluvial


def get_sensitivity_data(scenario_list, fname, SSPs, RCP="2p6", source="masked", only_first_scenario=False):
    """
    method to get data from different csvs scattered across different scenario definitions. Useful for collating results
    across the sensitivity analyses
    :param scenario_list: a list of all scenarios to be considered in the sensitivity analysis
    :param RCP: the RCP of the scenario being used. defaults to 2.6
    :param fname: the filename of the desired data sheet
    :param source: whether to get the masked or original source data
    :param SSPs: a list of SSPs for which to extract data
    :param only_first_scenario: only grabs the first scenario in the list for the sensitivity analysis
    :return:
    """
    if SSPs is None:
        SSPs = c.GCAMConstants.SSPs
    all_data = pd.DataFrame(columns=c.GCAMConstants.column_order)
    if only_first_scenario:
        scenario_list = [scenario_list[0]]
    for nonBaselineScenario in scenario_list:
        pyrolysis_df = pd.read_csv("data/gcam_out/" + str(nonBaselineScenario) + "/" + RCP + "/" + source + "/" + fname + ".csv")
        pyrolysis_df["Version"] = nonBaselineScenario
        pyrolysis_df = pyrolysis_df[pyrolysis_df[['SSP']].isin(SSPs).any(axis=1)]
        if all_data.empty:
            all_data = pyrolysis_df
        else:
            all_data = pd.concat([all_data, pyrolysis_df])
    return all_data


def drop_missing(df, max_length=100):
    """
    drops columns that contain "missing" from the dataframe for cleaner output .csv files
    :param df: dataframe with columns to be removed
    :return:
    """
    columns_to_drop = []
    for col in df.columns:
        if df[col].dtype == 'object':  # Check if the column is of type string (object)
            if any(df[col].astype(str).str.len() > max_length):
                columns_to_drop.append(col)

    df = df.drop(columns=columns_to_drop, axis=1)
    return df


def seq_C(row, product_column, modification_column):
    """
    calculates the sequestered C from the net change in C due to biochar
    :param row: the row of the dataframe being changed
    :param product_column: the column containing the identifying product
    :param modification_column: the column with data to be changed
    :return: the amount of net C that is sequestered
    """
    if row[product_column] in ["pork biochar", "goat biochar"]:
        return row[modification_column] * .5853 # from the ncomms spreadsheet
    elif row[product_column] in ["beef biochar", "dairy biochar"]:
        return row[modification_column] * .5391 # from the ncomms spreadsheet
    elif row[product_column] in ["poultry biochar"]:
        return row[modification_column] * .5316 # from the ncomms spreadsheet
    else:
        return 0


def avd_C(row, product_column, modification_column):
    """
    calculates the avoided C from the net change in C due to biochar
    :param row: the row of the dataframe being changed
    :param product_column: the column containing the identifying product
    :param modification_column: the column with data to be changed
    :return: the amount of net C that is sequestered
    """
    if row[product_column] in ["pork biochar", "goat biochar"]:
        return row[modification_column] * (1-.5853) # from the ncomms spreadsheet
    elif row[product_column] in ["beef biochar", "dairy biochar"]:
        return row[modification_column] * (1-.5391) # from the ncomms spreadsheet
    elif row[product_column] in ["poultry biochar"]:
        return row[modification_column] * (1-.5316) # from the ncomms spreadsheet
    else:
        return 0


def ghg_ER(row, product_column, modification_column):
    """
    calculates the avoided C from the net change in C due to biochar
    :param row: the row of the dataframe being changed
    :param product_column: the column containing the identifying product
    :param modification_column: the column with data to be changed
    :return: the amount of net C that is sequestered
    """
    # these values are already negative
    if row[product_column] in ["CH4"]:
        return row[modification_column] * 23  # emissions reduction, GWP from Ncomms spreadsheet, as all other ghg emissions reduction/CDR are negative, add a - sign to the returned value
    elif row[product_column] in ["N2O"]:
        return row[modification_column] * 296 # from the ncomms spreadsheet
    else:
        return 0



def avd_soil_emissions(row, product_column, modification_column):
    """
    calculates the avoided soil N2O emissions from biochar
    :param row: the row of the dataframe being changed
    :param product_column: the column containing the identifying product
    :param modification_column: the column with data to be changed
    :return: the amount of net C that is sequestered
    """
    if row['Version'] == "HighBiocharSoilN2O":
        temp = row[
                   modification_column] / .54  # counterfactual for full GHG emissions (parameter from ncomms spreadsheet)
        temp = temp * (
                1 - .54) * -296  # emissions reduction, GWP from Ncomms spreadsheet, as all other ghg emissions reduction/CDR are negative, add a - sign to the returned value
        return temp
    elif row['Version'] == "LowBiocharSoilN2O":
        temp = row[
                   modification_column] / 1.39  # counterfactual for full GHG emissions (parameter from ncomms spreadsheet)
        temp = temp * (
                1 - 1.39) * -296  # emissions reduction, GWP from Ncomms spreadsheet, as all other ghg emissions reduction/CDR are negative, add a - sign to the returned value
        return temp
    else:
        temp = row[
                   modification_column] / .98  # counterfactual for full GHG emissions (parameter from ncomms spreadsheet)
        temp = temp * (
                    1 - .98) * -296  # emissions reduction, GWP from Ncomms spreadsheet, as all other ghg emissions reduction/CDR are negative, add a - sign to the returned value
        return temp


def relabel_feeds(row):
    """
    relabels animal feeds for human clarity
    :param row: row in pandas dataframe
    :return: relabeled animal feed
    """
    if row['product'] in ["FodderHerb_Residue"]:
        return 'Fodder - Herb'
    elif row['product'] in ["FeedCrops"]:
        return 'Feed Crops'
    elif row['product'] in ["Pasture_FodderGrass"]:
        return 'Fodder - Grass'
    elif row['product'] in ["Scavenging_Other"]:
        return 'Scavenging'
    else:
        return "error"


def relabel_animals(row):
    """
    relabels animal products in GCAM for human clarity
    :param row: row in pandas dataframe
    :return: relabeled output
    """
    if row['product'] in ["SheepGoat"]:
        return "Sheep & Goat"
    else:
        return row['product']


def relabel_staple(row, column):
    """
        lambda function to relabel GCAM food categories for greater accessibility. From: A_demand_technology.csv
        :param row: row of data
        :param column: column of data to be processed
        :return: updated name of GCAM region
        """
    food = row[column]
    if food == "Corn" or food == "Wheat" or food == "OtherGrain" or food == "Rice" or food == "RootTuber":
        return "Staples"
    else:
        return "Non-Staples"


def get_CI(dataframe, products, alpha=0.95):
    """
    produce a confidence interval (defaults to 95%), returning the mean, median, and mode for each unique product in a column
    :param dataframe: dataframe with all the necessary data
    :param products: column of pandas dataf
    :return: a dataframe containing the min, median, and max for each product for all valid years
    """
    unique_products = dataframe[products].unique()
    lmu = pd.DataFrame(columns=dataframe.columns)
    for i in unique_products:
        # get a dataframe just for this product
        data = dataframe[dataframe[products] == i].copy(deep=True)

        # copy 3 rows over to lmu
        output_vals = data.head(3).copy(deep=True)

        if len(output_vals) == 3:
            output_vals["Version"] = ["Lower CI", "Median", "Upper CI"]

            # solve the CI for each year
            for j in c.GCAMConstants.future_x:
                np_data = data[str(j)].dropna().values  # get data for a particular year
                sMu = np.mean(np_data)
                median = np.median(np_data)
                sem = stats.sem(np_data)
                n = len(np_data)
                df = n - 1
                lower, upper = stats.t.interval(alpha, df=df, loc=sMu, scale=sem)  # confidence interval with equal areas around the mean

                # add lower, mean, upper to output dataframe
                output_vals[str(j)] = [lower, median, upper]

            # add lower level low mean upper dataframe to higher level one
            lmu = pd.concat([lmu, output_vals])

    # add the returned dataframe to the original frame as appended rows
    return lmu