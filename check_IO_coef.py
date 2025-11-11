import math
import sys
import numpy as np
import constants as c
import pandas as pd


def IO_check(comb, coefficients, assert_str, nonBaselineScenario, products):
    """
    Does a check over years and SSPs for
    :param comb: the dataframe containing data
    :param coefficients: the constant coefficients being analyzed
    :param assert_str: the string to be used in the error message
    :param nonBaselineScenario: the scenario being compared to the baseline scenario
    :param products: the product being analyzed
    :return: a mask of years and data inputs
    """
    to_return = list()
    for i in c.GCAMConstants.future_x:
        comb["check_" + str(i)] = comb[str(i) + "_x"] / comb[str(i) + "_y"]
    for k in c.GCAMConstants.SSPs:
        # check coef on an SSP basis
        comb_SSP = comb[comb[['SSP']].isin([k]).any(axis=1)]
        if len(comb_SSP) == 0:
            to_print = assert_str + " fails for product " + str(products) + " in " + str(
                k) + " because no data is available"
            print(to_print)
            for b in c.GCAMConstants.future_x:
                to_return.extend(tuple([[k, b]]))
        else:
            for i in c.GCAMConstants.future_x:
                try:  # if mean coefficient isn't within 5% of target
                    if not math.isinf(comb_SSP["check_" + str(i)].mean()) or not np.isnan(comb_SSP["check_" + str(i)].mean()):
                        assert (abs(coefficients[str(products)]) * .95 <
                                abs(comb_SSP["check_" + str(i)].mean()) <
                                abs(coefficients[str(products)]) * 1.05)
                    else:
                        print(assert_str + " is ignored for product " + str(products) + " in year " + str(
                            i) + " in " + str(
                            k) + \
                              " with mean " + str(comb_SSP["check_" + str(i)].mean()))
                except AssertionError:  # print out that the scenario is no good
                    to_print = assert_str + " fails for product " + str(products) + " in year " + str(i) + " in " + str(
                        k) + \
                               " with mean " + str(comb_SSP["check_" + str(i)].mean())
                    print(to_print)
                    to_return.extend(tuple([[k, i]]))
    return to_return


def getMask(nonBaselineScenario, RCP, filepath):
    """
    get the mask for the data pre-processing
    :param nonBaselineScenario: the non-baseline GCAM scenario being analyzed
    :param RCP: the RCP pathway for the non-baseline scenario
    :param filepath: filepath for log file of errors
    :return: the mask (or nothing if it is a released GCAM model) to update the datasets
    """
    old_stdout = sys.stdout
    with open(filepath + "mask_log.txt", "w") as log_file:
        sys.stdout = log_file
        mask = list()

        # read in data
        supply = pd.read_csv("data/gcam_out/" + str(
            nonBaselineScenario) + "/" + RCP + "/original/supply_of_all_markets.csv")  # current wd is /xml for some reason
        co2_emissions = pd.read_csv("data/gcam_out/" + str(
            nonBaselineScenario) + "/" + RCP + "/original/CO2_emissions_by_tech_excluding_resource_production.csv")

        for j in ["beef", "dairy", "goat", "pork", "poultry"]:
            # carbon yields from manure is the test parameter to check that the GCAM model doesn't have an error in its output
            products = [str(j) + " manure"]
            co2 = co2_emissions[co2_emissions[['technology']].isin([str(j) + "_biochar"]).any(axis=1)]
            manure = supply[supply['product'].str.contains("|".join(products))]

            # get global manure supply
            manure = manure.groupby(['SSP']).sum()
            manure['GCAM'] = 'Global'
            comb = pd.merge(co2, manure, how="inner", on=['GCAM', 'SSP'])
            if "LowCarbonStability" in filepath:
                mask.extend(x for x in
                            IO_check(comb, c.GCAMConstants.low_manure_C_ratio, "assert C and manure", nonBaselineScenario,
                                     str(j)) if x not in mask)
            elif "HighCarbonStability" in filepath:
                mask.extend(x for x in
                            IO_check(comb, c.GCAMConstants.high_manure_C_ratio, "assert C and manure", nonBaselineScenario,
                                     str(j)) if x not in mask)
            else:
                mask.extend(x for x in
                            IO_check(comb, c.GCAMConstants.manure_C_ratio, "assert C and manure", nonBaselineScenario,
                                     str(j)) if x not in mask)

        print(mask)
        sys.stdout = old_stdout
        return mask


if __name__ == '__main__':
    # edit the below lines if running the main program
    nonBaselineScenario = "test"
    RCP = "6p0"
    filepath = "data/gcam_out/" + nonBaselineScenario + "/" + RCP + "/"
    getMask(nonBaselineScenario, RCP, filepath)
