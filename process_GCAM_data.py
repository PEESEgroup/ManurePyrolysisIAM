import numpy as np
import pandas as pd
import csv
import constants as c
from itertools import islice
import check_IO_coef
import os


def split_file(fname):
    """
    splits a single GCAM csv file into seperate csv files, 1 for each query
    :param fname: the filename of the .csv file being processed
    :return: a dictionary of dataframes, 1 for each query, with the title as the key
    """
    # open file a first time to find the locations of the headers
    title_rows = []
    frames = {}
    names = {}
    row_counter = 0
    frames_key = ""

    with open(fname) as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            if len(row) == 1:  # all header rows have a length of 1
                if "java.lang.Exception" not in row[0]:
                    title_rows.append(row_counter)
            row_counter = row_counter + 1

    print(title_rows)

    # open the file a second time, this time reading chunks of data into dataframes based on locations identified above
    with open(fname) as f:
        for i in title_rows:
            if title_rows.index(i) + 1 == len(title_rows):
                # if we reached the last chunk in the list, just grab whatever remains
                line_chunk = [line for line in f]
            else:
                next_title = int(title_rows[title_rows.index(i) + 1])
                line_chunk = list(islice(f, next_title - i))

            # get title and column information from the .csv file
            if len(line_chunk) < 2:  # sometimes the model has an oopsie and doesn't have carbon sequestration data
                df = pd.DataFrame(columns=c.GCAMConstants.column_order)
            else:
                key = line_chunk[0].rstrip()
                col = line_chunk[1].split(",")
                # insert a new column into the list because when the data_list is broken into strings,
                # it breaks the scenario column into 2
                col.insert(1, "date")
                col[-1] = "newline_lmao"  # renaming last column because \n is a terrible column name

                # processing the data
                data_strings = line_chunk[2:]
                data_list = [row.split(",") for row in data_strings]

                # add data to dataframe
                if "depth" in data_list[0][5]:
                    col.insert(5, "depth")
                df = pd.DataFrame(data_list, columns=col)
                df = df.dropna()

                # add dataframe to dictionary
                # if the key already exists, ensure that they key is unique
                if key in names.keys():
                    frames_key = key + str(names[key])
                    names[key] = names[key] + 1
                else:
                    names[key] = 1
                    frames_key = key + "0"
            frames[frames_key] = df

    # combine files with multiple parts
    for item in names.items():
        # if there are more than one entry for that file
        append_df = pd.DataFrame()
        for i in range(item[1]):
            # add the GCAM region to the dataframe
            df = frames[item[0] + str(i)]
            df = df.loc[:, ~df.columns.duplicated()]

            # add the GCAM region to the master datadframe
            append_df = pd.concat([append_df, df], ignore_index=True)

            # remove the old dataframe from the dictionary
            del frames[item[0] + str(i)]

        frames[item[0]] = append_df
    return frames


def process_file(value, fname):
    """
    Ensures that each query .csv file has the same rows and columns - renames columns to be easier to graph later
    :param value: a dataframe containing the results of a query
    :param fname: the filepath to the location of the query, which stores version information
    :return: a standardized dataframe
    """
    # get model version
    words = fname.split("/")
    value['Version'] = words[2] + " " + words[3]

    # get GCAM market region information
    if 'market' in value.columns:
        value['GCAM'] = value.apply(lambda row: label_market_as_region(row), axis=1)
        value['product'] = value.apply(lambda row: label_market_as_product(row), axis=1)
    else:
        if 'region' in value.columns:
            value = value.rename(columns={'region': 'GCAM'})
        else:
            value['GCAM'] = c.GCAMConstants.missing
        value['product'] = c.GCAMConstants.missing

    # get SSP information
    if 'scenario' in value.columns:
        value['SSP'] = value.apply(lambda row: label_ssp(row), axis=1)
    else:
        value['SSP'] = c.GCAMConstants.missing

    # get sector, subsector, technology information, and make missing where possible
    if 'sector' not in value.columns:
        value['sector'] = c.GCAMConstants.missing
    if 'subsector' not in value.columns:
        value['subsector'] = c.GCAMConstants.missing
    if 'technology' not in value.columns:
        value['technology'] = c.GCAMConstants.missing
    if 'output' not in value.columns:
        value['output'] = c.GCAMConstants.missing
    if 'concentration' not in value.columns:
        value['concentration'] = c.GCAMConstants.missing
    if 'input' not in value.columns:
        value['input'] = c.GCAMConstants.missing
    if "LandLeaf" not in value.columns:
        value['LandLeaf'] = c.GCAMConstants.missing
    if "GHG" not in value.columns:
        value['GHG'] = c.GCAMConstants.missing
    if "fuel" not in value.columns:
        value['fuel'] = c.GCAMConstants.missing

    # reorder columns
    return value[c.GCAMConstants.column_order]


def label_ssp(row):
    """
    labels each row with SSP informtion extracted from the scenario column
    :param row: a row in a dataframe
    :return: the SSP number as "SSP#", or a default value if it is not a valid SSP number
    """
    for i in c.GCAMConstants.SSPs:
        if i in row['scenario']:
            return i
    return c.GCAMConstants.missing


def label_market_as_region(row):
    """
    labels each row with regional information extracted from the market column
    :param row: a row in the datafame
    :return: the GCAM region name, or a default value if it is not a region name
    """
    for j in c.GCAMConstants.GCAM_region:
        if j in row['market']:
            return j
    return c.GCAMConstants.missing


def label_market_as_product(row):
    """
    labels the key products we are interested in as extracted from the market column
    :param row: a row in the dataframe
    :return: the key market, or a default value if it is not a key market
    """
    for j in c.GCAMConstants.GCAM_region:
        if j in row['market']:
            return row['market'].replace(j, '')
    return c.GCAMConstants.missing


def masking(dataframe, mask):
    """
    masks errors in input data with np.nan
    :param dataframe: dataframe to be masked
    :param mask: list of SSP-year pairs with model errors
    :return: dataframe with relabled rows
    """
    for i in mask:
        year = str(i[1])
        SSP = str(i[0])
        # if SSP is not in the dataframe, then there's no need to do the computations to apply the mask
        length_SSP = dataframe[dataframe[['SSP']].isin([SSP]).any(axis=1)]
        if len(length_SSP) > 0:
            dataframe.loc[:, str(i[1])] = dataframe.apply(lambda row: apply_mask(row, year, SSP), axis=1)
    return dataframe


def apply_mask(row, year, SSP):
    """
    For a given row, apply a mask if SSP in row matches SSP with error
    :param row: a row from a pd dataframe
    :param year: the year to be masked
    :param SSP: the SSP to be masked
    :return: np.nan if the SSP matches, otherwise preserve original value
    """
    if row["SSP"] == SSP:
        return np.nan
    else:
        return row[year]


def main():
    """
    control block for this file
    :return: nothing, but writes out .csv files to a relative directory
    """
    for i in c.GCAMConstants.GCAMDB_filenames:  # first time through process the original data
        csvs = split_file(i)  # split file based on header rows

        # create directories if they don't already exist
        dir_path = i.split("/")  # fix the filename
        dir_path[-1] = "original/"
        original_fname = "/".join(dir_path)
        if not os.path.exists(original_fname):
            os.makedirs(original_fname)

        for item in csvs.items():  # for each file
            df = process_file(item[1], i)  # preprocess the data
            dir_path = i.split("/")  # fix the filename
            dir_path[-1] = "original/" + str(item[0]) + ".csv"
            new_fname = "/".join(dir_path)
            new_fname = new_fname.replace(")", "").replace("(", "").replace("\\", "").replace(" ", "_").replace("b/t",
                                                                                                                "between")
            print(new_fname)
            df.to_csv(new_fname, index=False)  # save the original file

    for i in c.GCAMConstants.GCAMDB_filenames:  # second time through process the mask
        csvs = split_file(i)  # split file based on header rows

        # get the mask if it isn't the reference data
        if i.split("/")[2] != "ref" and i.split("/")[2] != "released":
            mask = check_IO_coef.getMask(i.split("/")[2], i.split("/")[3], "/".join(i.split("/")[:-1])+"/")

            # create directories if they don't already exist
            dir_path = i.split("/")  # fix the filename
            dir_path[-1] = "masked/"
            masked_fname = "/".join(dir_path)
            if not os.path.exists(masked_fname):
                os.makedirs(masked_fname)
                print("made masked directory")

            for item in csvs.items():  # for each file
                df = process_file(item[1], i)  # preprocess the data
                df = masking(df, mask)  # apply the mask to the dataframe
                dir_path = i.split("/")  # fix the filename
                dir_path[-1] = "masked/" + str(item[0]) + ".csv"
                new_fname = "/".join(dir_path)
                new_fname = (new_fname.replace(")", "").replace("(", "").replace("\\", "").
                             replace(" ", "_").replace("b/t", "between"))
                print(new_fname)
                df.to_csv(new_fname, index=False)  # save the masked file


if __name__ == '__main__':
    main()
