import geopandas as gpd
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import pandas as pd
import os
import constants as c
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pylab import *
import data_manipulation
import plotly.express as px


def world_data(data):
    """
    merges the data to be plotted into a geopandas dataframe and selects relevant columns for plotting
    :param data: the dataframe containing GCAM data to be plotted
    :return: a pandas dataframe containing relevant year and column data
    """
    world = gpd.read_file(c.GCAMConstants.processed_map_loc)
    merged = pd.merge(world, data, on="GCAM", how='left')
    merged = merged.replace(c.GCAMConstants.missing, np.nan)

    return merged


def basin_data(data, column, title):
    """
    merges the data to be plotted into a geopandas dataframe and selects relevant columns for plotting
    :param data: the dataframe containing GCAM data to be plotted
    :return: a pandas dataframe containing relevant year and column data
    """
    # read in data
    data["GLU"] = data["GLU"].str.replace("GLU00", "").str.replace("GLU0", "").str.replace("GLU", "") # remove GLU codes and leading 0s
    data["GLU"] = data["GLU"].astype("int64")
    data["SSP"] = "SSP1"
    basins = gpd.read_file(c.GCAMConstants.basin_map_loc)
    n_products = len(data["GCAM_subsector"].unique())

    # set up plot info
    cmap = plt.colormaps.get_cmap('viridis')
    normalizer = Normalize(min(data[column]), max(data[column]))
    im = cm.ScalarMappable(norm=normalizer, cmap=cmap)

    # for each crop in the subsector
    for i in data["GCAM_subsector"].unique():
        fig, axs = plt.subplots(1, 1, sharex='all', sharey='all', gridspec_kw={'wspace': 0.2, 'hspace': 0.2})
        crop = data[data["GCAM_subsector"] == str(i)]
        merged = pd.merge(basins, crop, left_on=["glu_id", "reg_id"], right_on=["GLU","GCAM_region_ID"], how='left')
        merged = merged.replace(c.GCAMConstants.missing, np.nan)

        if str(i) == "SugarCropC4":
            subplot_title = title + " "  + "Sugar Cane"
        elif str(i) == "SugarCrop":
            subplot_title = title + " " + "Sugar Beet and Sugar Crops"
        else:
            subplot_title = title + " "  + str(i)

        plot_world_on_axs(
            map_plot=merged,
            axs=axs,
            cmap=cmap,
            counter=0,
            plot_title=subplot_title,
            plotting_column=column,
            ncol=1,
            nrow=1,
            normalizer=normalizer)
        units = "kg/ha"

        # update the figure with shared colorbar
        lab = units
        cax = fig.add_axes([0.08, 0.20, 0.02, 0.45])
        fig.colorbar(im, cax=cax, shrink=0.5, orientation="vertical", label=lab)

        # change figure size and dpi
        fig.set_dpi(300)
        plt.savefig("data/data_analysis/images/maps/" + subplot_title + ".png", dpi=300)
        plt.show()
        merged.drop("geometry",axis=1).to_csv("data/data_analysis/supplementary_tables/maps/" + subplot_title + ".csv")


def get_subplot_dimensions(list_products):
    """
    returns the number of subplot dimensions based on the number of items being plotted
    :param list_products: the list of things being plotted on one plot
    :return: the number of rows and the number of columns for the plot
    """
    if len(list_products) == 1:
        return 1, 1
    elif len(list_products) == 2:
        return 2, 1
    elif len(list_products) == 3:
        return 2, 2
    elif len(list_products) == 4:
        return 2, 2
    elif len(list_products) == 5 or len(list_products) == 6:
        return 2, 3
    elif len(list_products) == 7 or len(list_products) == 8:
        return 2, 4
    elif len(list_products) == 9:
        return 3, 3
    elif len(list_products) in [10, 11, 12]:
        return 3, 4
    elif len(list_products) in [13, 14, 15, 16]:
        return 4, 4
    elif len(list_products) in [17, 18, 19, 20]:
        return 4, 5
    else:
        raise ValueError("too many products. Can only plot 20 products at a time")


def plot_world(dataframe, products, SSPS, groupby, column, years, title, RCP, nonBaselineScenario):
    """
    control function for plotting any plot that is placed on a world map
    :param years: the years for which the plot is to be evaluated
    :param dataframe: the dataframe containing data to be plotted
    :param products: the list of products that are to be plotted
    :param SSPS: the list of SSPs by which the data is to be plotted
    :param groupby: how the data should be grouped. Accepted values include "SSP", "product", "year"
    :param column: the column on which the data should be filtered by product
    :param title: the title of the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: shows the relevant plot
    """
    if groupby == "SSP":
        plot_world_by_SSP(dataframe, products, column, years, SSPS, title, RCP, nonBaselineScenario)
    elif groupby == "product":
        plot_world_by_products(dataframe, products, column, years, SSPS, title, RCP, nonBaselineScenario)
    elif groupby == "year":
        plot_world_by_years(dataframe, products, column, years, SSPS, title, RCP, nonBaselineScenario)
    else:
        raise ValueError("only 'SSP', 'product', and 'year' are considered valid groupings at this time")


def plot_world_by_SSP(dataframe, products, column, year, SSP, title, RCP, nonBaselineScenario):
    """
    For a given product, plot its values in all SSPs for a given year
    :param dataframe: the dataframe containing the data to be plotted
    :param products: the product being plotted
    :param column: the column in the dataframe containing the product
    :param year: the year being chosen
    :param title: the title for the plot
    :param SSP: the set of shared socioeconomic pathways being used as evaluation
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    """
    for j in year:
        for i in products:
            try:
                counter = 0
                units = "N/A"
                # get plot information
                axs, cmap, fig, im, ncol, normalizer, nrow = create_subplots(
                    dataframe=dataframe,
                    inner_loop_set=SSP,
                    products=[i],
                    year=[j],
                    SSP=SSP,
                    product_column=column,
                    title=title)

                for k in SSP:
                    subplot_title = str(k)
                    units = get_df_to_plot(
                        dataframe=dataframe,
                        ncol=ncol,
                        nrow=nrow,
                        fig=fig,
                        axs=axs,
                        cmap=cmap,
                        normalizer=normalizer,
                        counter=counter,
                        column=column,
                        products=i,
                        SSPs=k,
                        years=j,
                        subplot_title=subplot_title)
                    counter = counter + 1

                # update the figure with shared colorbar
                dl = len(SSP)
                lab = units
                add_colorbar_and_plot(axs, dl, fig, im, lab, ncol, nrow, title, RCP, nonBaselineScenario)
            except ValueError as e:
                print(e)


def plot_world_by_products(dataframe, products, column, year, SSP, title, RCP, nonBaselineScenario):
    """
    For each SSP, plots all relevant products
    :param SSP: the SSP scenario for the plot
    :param year: The year of data to be plotted
    :param dataframe: dataframe containing the data to be plotted
    :param products: the data to be plotted
    :param column: the column for which the data is to be filtered
    :param title: the title for the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: shows the relevant plot
    """
    for j in year:
        for k in SSP:
            try:
                counter = 0
                units = "N/A"
                # get plot information
                axs, cmap, fig, im, ncol, normalizer, nrow = create_subplots(
                    dataframe=dataframe,
                    inner_loop_set=products,
                    products=products,
                    year=[j],
                    SSP=[k],
                    product_column=column,
                    title=title)

                # iterate through all subplots
                for i in products:
                    subplot_title = str(i)
                    units = get_df_to_plot(
                        dataframe=dataframe,
                        ncol=ncol,
                        nrow=nrow,
                        fig=fig,
                        axs=axs,
                        cmap=cmap,
                        normalizer=normalizer,
                        counter=counter,
                        column=column,
                        products=i,
                        SSPs=k,
                        years=j,
                        subplot_title=subplot_title)
                    counter = counter + 1

                # update the figure with shared colorbar
                dl = len(products)
                lab = units
                add_colorbar_and_plot(axs, dl, fig, im, lab, ncol, nrow, title, RCP, nonBaselineScenario)
            except ValueError as e:
                print(e)


def axs_params(ax, plot_title):
    """
    Provides uniform formatting for all axes
    :param ax: the axis being formatted
    :param plot_title: the title of the subplot
    :return: N/A
    """
    ax.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    ax.margins(x=0.005, y=0.005)
    ax.set_title(str(plot_title), fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


def plot_world_by_years(dataframe, products, column, year, SSP, title, RCP, nonBaselineScenario):
    """
    For each SSP, plots all relevant products
    :param SSP: the SSP scenario for the plot
    :param year: The year of data to be plotted
    :param dataframe: dataframe containing the data to be plotted
    :param products: the data to be plotted
    :param column: the column for which the data is to be filtered
    :param title: the title for the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: shows the relevant plot
    """
    for k in SSP:
        for i in products:
            try:
                counter = 0
                units = "N/A"
                # get plot information
                axs, cmap, fig, im, ncol, normalizer, nrow = create_subplots(
                    dataframe=dataframe,
                    inner_loop_set=year,
                    products=[i],
                    year=year,
                    SSP=[k],
                    product_column=column,
                    title=title)

                # iterate through all subplots
                for j in year:
                    subplot_title = str(j)
                    units = get_df_to_plot(
                        dataframe=dataframe,
                        ncol=ncol,
                        nrow=nrow,
                        fig=fig,
                        axs=axs,
                        cmap=cmap,
                        normalizer=normalizer,
                        counter=counter,
                        column=column,
                        products=i,
                        SSPs=k,
                        years=j,
                        subplot_title=subplot_title)
                    counter = counter + 1

                # update the figure with shared colorbar
                dl = len(year)
                lab = str(units)
                add_colorbar_and_plot(axs, dl, fig, im, lab, ncol, nrow, title, RCP, nonBaselineScenario)
            except ValueError as e:
                print(e)


def get_df_to_plot(dataframe, ncol, nrow, fig, axs, cmap, normalizer, counter, column, products, SSPs, years,
                   subplot_title):
    """
    This method filters data and prepare it for plotting
    :param dataframe: the data being plotted
    :param ncol: the number of columns of subplots
    :param nrow: the number of rows of subplots
    :param fig: matplotlib fig
    :param axs: the matplotlib axs being plotted
    :param cmap: the colormap being used
    :param normalizer: the normalizer bieng used
    :param counter: the plot counter
    :param column: the column on which the products are located
    :param products: the products being plotted
    :param SSPs: the SSPs being plotted
    :param years: the years being plotted
    :param subplot_title: the title for the subplot
    :return: unit label for the subplot
    """
    if column != "":
        filter_data = dataframe[dataframe[[column]].isin([products]).any(axis=1)]
        filter_data = filter_data[filter_data[['SSP']].isin([SSPs]).any(axis=1)]
    else:
        filter_data = dataframe

    # if there is no data in the filter data, delete all following axis
    if filter_data.empty:
        if ncol == 1:
            if nrow == 1:
                fig.delaxes(axs)
            else:
                fig.delaxes(axs[counter])
        else:
            fig.delaxes(axs[int(counter / ncol), int(counter % ncol)])
    else:
        map_to_plot = world_data(filter_data)
        plot_world_on_axs(
            map_plot=map_to_plot,
            axs=axs,
            cmap=cmap,
            counter=counter,
            plot_title=subplot_title,
            plotting_column=years,
            ncol=ncol,
            nrow=nrow,
            normalizer=normalizer)
        units = map_to_plot['Units'].unique()
        for i in units:
            if str(i) != "nan":
                return str(i)


def plot_world_on_axs(map_plot, axs, cmap, counter, plot_title, plotting_column, ncol, normalizer, nrow):
    """
    Plot the world map in a subplot
    :param map_plot: the dataframe containing map info and the data being plotted
    :param axs: the subplot axis on which the map is being plotted
    :param cmap: the color map for the choropleth map
    :param counter: the counter for which subplot is being counted
    :param plot_title: the subplot title
    :param plotting_column: the column containing data that is the basis of the choropleth map
    :param ncol: the number of columns of subplots
    :param normalizer: the data normalizer for the colorbar
    :param nrow: the number of rows of subplots
    :return: N/A
    """
    if ncol == 1:
        if nrow == 1:
            map_plot.plot(
                column=str(plotting_column),
                missing_kwds=dict(color='grey', label='No Data'),
                ax=axs,
                cmap=cmap,
                norm=normalizer)
            axs_params(axs, plot_title)
        else:
            map_plot.plot(
                column=str(plotting_column),
                missing_kwds=dict(color='grey', label='No Data'),
                ax=axs[int(counter / ncol)],
                cmap=cmap,
                norm=normalizer)
            axs_params(axs[int(counter / ncol)], plot_title)
    else:
        map_plot.plot(
            column=str(plotting_column),
            missing_kwds=dict(color='grey', label='No Data'),
            ax=axs[int(counter / ncol), int(counter % ncol)],
            cmap=cmap,
            norm=normalizer)
        axs_params(axs[int(counter / ncol), int(counter % ncol)], plot_title)


def create_subplots(dataframe, inner_loop_set, products, year, SSP, product_column, title):
    """
    Creates the boilerplate subplots and colorbars
    :param inner_loop_set: the list of products being iterated over for the different subplots
    :param SSP: The SSP of the scenario
    :param product_column: The column on which the different products being graphed vary
    :param dataframe: the dataframe being evaluated
    :param products: the list of products being evaluated
    :param year: the year in which the data is plotted
    :param title: the title of the plot
    :return: a set of figure objects
    """
    # at this stage, if this df is empty, then we know that there is no material to plot
    if dataframe.empty:
        raise ValueError("These products" + str(products) + "do not exist in this dataframe")
    df = dataframe.replace([np.inf, -np.inf], np.nan)

    # get subplot size
    nrow, ncol = get_subplot_dimensions(inner_loop_set)

    # make plots and color scheme
    fig, axs = plt.subplots(nrow, ncol, sharex='all', sharey='all', gridspec_kw={'wspace': 0.2, 'hspace': 0.2})
    cmap = plt.colormaps.get_cmap('viridis')
    fig.suptitle(title)
    normalizer = Normalize(min(df[str(i)].min() for i in year), max(df[str(i)].max() for i in year))
    im = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    plt.xticks(range(min([int(y) for y in year]), max([int(y) for y in year])+1, 5))
    return axs, cmap, fig, im, ncol, normalizer, nrow


def add_colorbar_and_plot(axs, datalength, fig, im, lab, ncol, nrow, fname, RCP, nonBaselineScenario):
    """
    Adds the colorbar to the graph and produces output
    :param axs: matplotlib axes
    :param datalength: the length of the data being plotted
    :param fig: matplotlib fig
    :param im: matplotlib im
    :param lab: label for the colorbar
    :param ncol: number of rows of axes
    :param nrow: number of columns of axes
    :param fname: filename for output image
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: plotted figure
    """
    del_axes_counter = 0
    for i in range(nrow * ncol - datalength):
        fig.delaxes(axs[int((datalength + i) / ncol), int((datalength + i) % ncol)])
        del_axes_counter = del_axes_counter + 1
    if nrow * ncol > datalength:
        if del_axes_counter == 3:
            axins1 = inset_axes(
                axs[int((datalength + 1) / ncol), int((datalength + 1) % ncol)],
                width=str(100 * del_axes_counter) + "%",  # width: 50% of parent_bbox width
                height="10%",  # height: 5%
                loc="center",
            )
        else:
            axins1 = inset_axes(
                axs[int(datalength / ncol), int(datalength % ncol)],
                width=str(100 * del_axes_counter) + "%",  # width: 50% of parent_bbox width
                height="10%",  # height: 5%
                loc="center",
            )
        axins1.xaxis.set_ticks_position("bottom")
        axins1.set_title(lab)
        fig.colorbar(im, cax=axins1, orientation="horizontal")
    elif ncol == 1:
        if nrow == 1:
            cax = fig.add_axes([0.9, 0.25, 0.02, 0.5])
            fig.colorbar(im, cax=cax, shrink=0.6, orientation="vertical", label=lab)
        else:
            fig.colorbar(im, ax=axs[:], shrink=0.6, orientation="vertical", label=lab)
    else:
        fig.colorbar(im, ax=axs[nrow - 1, :], shrink=0.6, orientation="horizontal", label=lab)

    # change figure size and dpi
    fig.set_dpi(300)
    if nrow * ncol == 6:
        fig.set_size_inches(12, 4)
    elif nrow * ncol == 8:
        fig.set_size_inches(16, 5.1)
    elif nrow * ncol == 20:
        fig.set_size_inches(16, 9)
    plt.savefig("data/data_analysis/images/"  + str(RCP) + "/"  +  fname + ".png", dpi=300)
    plt.show()


def plot_line_by_SSP(dataframe, products, column, SSP, differentiator, title, RCP, nonBaselineScenario):
    """
    plots a line graph with different subplot for each SSP
    :param dataframe: the dataframe with data to be plotted
    :param products: the list of products to be plotted
    :param column: the column i nthe dataframe in which the products can be found
    :param SSP: the SSP scenario
    :param differentiator: the column containing unique model names
    :param title: the title of hte plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    try:
        # get plot information
        axs, cmap, fig, im, ncol, normalizer, nrow = create_subplots(
            dataframe=dataframe,
            inner_loop_set=products,
            products=products,
            year=c.GCAMConstants.biochar_x,
            SSP=SSP,
            product_column=column,
            title=title)

        # find the number of model versions
        # get color scheme based on number of model versions
        versions = dataframe[differentiator].unique()
        colors, num_colors = get_colors(len(versions))
        counter = 0
        for i in products:
            color_counter = 0
            for k in SSP:
                sub_color = 0
                y = dataframe[(dataframe[column] == i) & (dataframe['SSP'] == k)]

                if not y.empty:
                    # plot all versions in y
                    for j in y[differentiator].unique():
                        # get y label
                        if len(versions) > 1:
                            lab = str(k) + str(j)
                        else:
                            lab = str(k)

                        # get color
                        color = colors[color_counter * num_colors + sub_color]
                        sub_color = sub_color + 1

                        # get line of data to plot and plot it
                        df = y[y[differentiator] == j]
                        y_to_plot = df.values.tolist()[0][c.GCAMConstants.skip_years:c.GCAMConstants.skip_years + len(
                            c.GCAMConstants.biochar_x)]  # only take the x values
                        plot_line_on_axs(
                            x=c.GCAMConstants.biochar_x,
                            y=y_to_plot,
                            lab=lab,
                            color=color,
                            axs=axs,
                            nrow=nrow,
                            ncol=ncol,
                            counter=counter)

                    # get units
                    units = y['Units'].unique()[0]
                    color_counter = color_counter + 1

            l, h = finalize_line_subplot(axs, units, str(i), ncol, nrow, counter)

            counter = counter + 1
        finalize_line_plot(fig, h, l, axs, nrow, ncol, counter, title, RCP, nonBaselineScenario)

    except ValueError as e:
        print(e)


def finalize_line_plot(fig, handles, labels, axs, nrow, ncol, counter, title, RCP, nonBaselineScenario):
    """
    adds a legend to the plot and removes unnecessary axes
    :param fig: matplotlib figure
    :param handles: matplotlib handles
    :param labels: matplotlib labels
    :param axs: maxplotlib axs
    :param nrow: the number of rows of subplots
    :param ncol: the number of columns of subplots
    :param counter: the number of used subplots
    :param title: the title of the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    if nrow * ncol == 1:
        fig.legend(handles, labels, bbox_to_anchor=(1, .895), facecolor='white', framealpha=1)
        plt.subplots_adjust(bottom=0.4, right=.7)
    elif nrow * ncol == 2:
        fig.legend(handles, labels, bbox_to_anchor=(1, .895), facecolor='white', framealpha=1)
        plt.subplots_adjust(bottom=0.4, right=.7)
    elif nrow * ncol == 4:
        fig.legend(handles, labels, bbox_to_anchor=(0.6, 0.6), facecolor='white', framealpha=1)
        plt.tight_layout()
    elif nrow * ncol == 6 and counter == 5:
        fig.legend(handles, labels, bbox_to_anchor=(0.9, 0.4), facecolor='white', framealpha=1)
        plt.tight_layout()
    elif nrow * ncol == 6:
        fig.legend(handles, labels, bbox_to_anchor=(0.15, 0.54), facecolor='white', framealpha=1)
        plt.tight_layout()
    else:
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 1), facecolor='white', framealpha=1)
        plt.tight_layout()

    # remove unnecessary axes
    for i in range(nrow * ncol - counter):
        fig.delaxes(axs[int((counter + i) / ncol), int((counter + i) % ncol)])

    plt.savefig("data/data_analysis/images/"  + str(RCP) + "/"  +  title + ".png", dpi=300)
    plt.show()


def plot_line_on_axs(x, y, lab, color, axs, nrow, ncol, counter):
    """
    Plots a line on the axis
    :param x: x values for the plot
    :param y: y values for the plot
    :param lab: label for the line
    :param color: the color of the line
    :param axs: the axis being plotted on
    :param nrow: the number of rows of axes
    :param ncol: the number of columns of axes
    :param counter: the number of the current subplot
    :return: N/A
    """
    if ncol == 1:
        if nrow == 1:
            axs.plot(x, y, label=lab, color=color)
        else:
            axs[int(counter / ncol)].plot(x, y, label=lab, color=color)
    else:
        axs[int(counter / ncol), int(counter % ncol)].plot(x, y, label=lab, color=color)


def finalize_line_subplot(axs, ylabel, title, ncol, nrow, counter):
    """
    Applies formatting to a subplot
    :param axs: the axis of the subplot
    :param ylabel: the label for the y-axis
    :param title: the title of the graph
    :param ncol: the number of rows of subplots
    :param nrow: the number of columns of subplots
    :param counter: the current subplot
    :return: labels and handles of the subplot
    """
    if ncol == 1:
        if nrow == 1:
            axs.set_ylabel(ylabel)
            axs.set_xlabel("Year")
            axs.set_title(title)
            handles, labels = axs.get_legend_handles_labels()
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        else:
            axs[int(counter / ncol)].set_ylabel(ylabel)
            axs[int(counter / ncol)].set_xlabel("Year")
            axs[int(counter / ncol)].set_title(title)
            handles, labels = axs[0].get_legend_handles_labels()
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    else:
        axs[int(counter / ncol), int(counter % ncol)].set_ylabel(ylabel)
        axs[int(counter / ncol), int(counter % ncol)].set_title(title)
        axs[int(counter / ncol), int(counter % ncol)].set_xlabel("Year")
        handles, labels = axs[0, 0].get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))

    return labels, handles


def plot_line_by_product(dataframe, products, column, SSP, differentiator, title, RCP, nonBaselineScenario):
    """
    Plots a line grouped by product
    :param dataframe: the data being plotted
    :param products: the list of products in the data being plotted
    :param column: the column by which the products can be differentiated
    :param SSP: the SSPs being plotted
    :param differentiator: a secondary column to differentiate the products
    :param title: the title of the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    # get plot information
    axs, cmap, fig, im, ncol, normalizer, nrow = create_subplots(
        dataframe=dataframe,
        inner_loop_set=SSP,
        products=products,
        year=c.GCAMConstants.biochar_x,
        SSP=SSP,
        product_column=column,
        title=title)

    # find the number of model versions
    # get color scheme based on number of model versions
    versions = dataframe[differentiator].unique()
    colors, num_colors = get_colors(len(versions))
    counter = 0
    try:
        for k in versions:
            color_counter = 0
            for i in products:
                y = dataframe[(dataframe[column] == i) & (dataframe[differentiator] == k)]

                if not y.empty:
                    # plot all versions in y
                    color = colors[color_counter]

                    # get line of data to plot and plot it
                    lower = y.values.tolist()[0][c.GCAMConstants.skip_years:c.GCAMConstants.skip_years + len(
                        c.GCAMConstants.biochar_x)]
                    y_to_plot = y.values.tolist()[1][c.GCAMConstants.skip_years:c.GCAMConstants.skip_years + len(
                        c.GCAMConstants.biochar_x)]  # only take the x values
                    upper = y.values.tolist()[2][c.GCAMConstants.skip_years:c.GCAMConstants.skip_years + len(
                        c.GCAMConstants.biochar_x)]
                    plot_line_on_axs(c.GCAMConstants.biochar_x, y_to_plot, str(i), color, axs, nrow, ncol, counter)
                    axs.fill_between(c.GCAMConstants.biochar_x, lower, upper, color=color, alpha = 0.2)

                    # get units
                    units = y['Units'].unique()[0]
                    l, h = finalize_line_subplot(axs, units, str(k), ncol, nrow, counter)
                    color_counter = color_counter + 1

            counter = counter + 1
        finalize_line_plot(fig, h, l, axs, nrow, ncol, counter, title, RCP, nonBaselineScenario)

    except ValueError as e:
        print(e)


def get_colors(num_versions):
    """
    gets a color mapping based on the number of versions of the product
    :param num_versions: the number of unique entries for each product
    :return: a list containing the requisite colors, the number of colors for each product
    """
    if num_versions == 1:
        num_sub_colors = 1
        return ["#BFBE43", "#74A751", "#698FC6", "#DD9452", "#C16861", "#A577A8", "#72B1B4", "#DCC060", "#AD9077",
                "#9299A9"], num_sub_colors
    elif num_versions == 2:
        cmap = matplotlib.colormaps.get_cmap('tab20')
        num_sub_colors = 2
    elif num_versions == 4:
        cmap = matplotlib.colormaps.get_cmap('tab20b')
        num_sub_colors = 4
    else:
        return ["#1F78C8", "#ff0000", "#33a02c", "#6A33C2", "#ff7f00", "#565656",
                "#FFD700", "#a6cee3", "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F",
                "#999999", "#EEE685", "#C8308C", "#FF83FA", "#C814FA", "#0000FF",
                "#36648B", "#00E2E5", "#00FF00", "#778B00", "#BEBE00", "#8B3B00",
                "#A52A3C"], 25
    return [matplotlib.colors.rgb2hex(c) for c in cmap.colors], num_sub_colors


def plot_stacked_bar_product(df, year, SSP, column, title, RCP, nonBaselineScenario):
    """
    Plots a stacked bar graph
    :param df: dataframe
    :param year: list of year to be plotted
    :param SSP: list of SSPs to be plotted
    :param column: column in the df to be plotted
    :param title: title for the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    try:
        # get subplot information
        nrow, ncol = get_subplot_dimensions([year])
        fig, axs = plt.subplots(nrow, ncol, sharex='all', sharey='all', gridspec_kw={'wspace': 0.2, 'hspace': 0.2})
        colors, num_colors = get_colors(1)

        # format table
        plot_df = df[df[['SSP']].isin(SSP).any(axis=1)]
        if not isinstance(year, list):
            plot_df = plot_df.loc[:, [year, 'SSP', 'GCAM', column]]
            plot_df = plot_df.pivot(index='GCAM', columns=column, values=year)
            plot_df = plot_df.loc[:, (plot_df != 0).any(axis=0)]
            # add a column to df for sorting, then remove it
            plot_df["sum"] = plot_df.abs().sum(axis=1)
            plot_df = plot_df.sort_values(by="sum", ascending=False)
            plot_df = plot_df.drop(['sum'], axis=1)
            # plot stacked bar chart
            plot_df.plot(kind="bar", stacked=True, color=colors, ax=axs)
        else:
            plot_df = pd.melt(plot_df, id_vars=['LandLeaf'], value_vars=[str(j) for j in year])
            plot_df = plot_df.pivot(index="variable", columns="LandLeaf", values="value")
            # plot stacked bar chart
            plot_df.plot(kind="bar", stacked=True, color=colors, ax=axs)

        # format plot
        axs.set_title(title)
        axs.set_ylabel(df["Units"].unique()[0])
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(bottom=0.5, right=.7, left=.15)
        plt.xticks(rotation=60, ha='right')
        ymin, ymax = axs.get_ylim()
        axs.set_ylim(-ymax, ymax)
        if title == "global land use change by year":
            plt.gcf().set_size_inches(7, 8)
        else:
            plt.gcf().set_size_inches(12, 8)
        plt.savefig("data/data_analysis/images/" + str(RCP) + "/"  + title + ".png", dpi=300)
        plt.show()

    except ValueError as e:
        print(e)


def plot_regional_vertical(dataframe, year, SSPs, y_label, title, x_column, y_column, x_label, RCP, nonBaselineScenario):
    """
    Plots regional data in a categorical scatterplot
    :param dataframe: data being plotted
    :param year: year of y_axis data
    :param SSPs: SSPs being evaluated
    :param y_label: ylabel for graph
    :param x_label: x-label for graph
    :param title: title of graph
    :param x_column: column used on the x-axis (formerly "GCAM")
    :param y_column: column used on the y-axis (formerly "columnn")
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    # get colors
    colors, divisions = get_colors(len(dataframe[y_column].unique()))

    # plot for each SSP
    for i in SSPs:
        dataframe = dataframe[dataframe['SSP'].str.contains(i)]
        for idx, item in enumerate(dataframe[y_column].unique()):
            df = dataframe.loc[dataframe[y_column] == str(item)]
            # scatter points
            plt.scatter(x=df[x_column], y=df[str(year)], color=colors[idx], label=str(item))

            # plot averages
            # plt.axhline(y=df[str(year)].mean(), color=colors[idx], linestyle='dashed', label=str(item) + " average")

        # finalize plot
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        plt.xticks(rotation=60, ha='right')
        plt.title(title)
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(bottom=0.4, right=.7)
        plt.show()
        plt.savefig("data/data_analysis/images/" + str(RCP) + "/"  +  title + ".png", dpi=300)


def plot_regional_vertical_avg(prices, year, SSPs, y_label, title, column, supply, RCP, nonBaselineScenario):
    """
    Plots regional data in a categorical scatterplot
    :param prices: price data being plotted
    :param year: evaluation column
    :param SSPs: SSPs being evaluated
    :param y_label: ylabel for graph
    :param title: title of graph
    :param column: column used to identify unique categories
    :param supply: supply data being plotted for weighted averages
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    # get colors
    colors, divisions = get_colors(5)

    # plot for each SSP
    printing_str = ""
    for i in SSPs:
        dataframe = prices[prices['SSP'].str.contains(i)]
        for idx, item in enumerate(dataframe[column].unique()):
            df_price = dataframe.loc[dataframe[column] == str(item)]

            # scatter points
            global_avg = df_price[df_price[['GCAM']].isin(["Global"]).any(axis=1)]
            df_price = df_price[~df_price[['GCAM']].isin(["Global"]).any(axis=1)]

            df_lower = df_price[df_price["Version"] == "Lower CI"]
            df_upper = df_price[df_price["Version"] == "Upper CI"]
            df_median = df_price[df_price["Version"] == "Median"]
            global_avg = global_avg[global_avg["Version"] == "Median"]

            plt.scatter(x=df_median["GCAM"], y=df_median[str(year)], edgecolors=colors[idx], facecolors='none', label=str(item))
            plt.scatter(x=df_lower["GCAM"], y=df_lower[str(year)], color=colors[idx], marker="*")
            plt.scatter(x=df_upper["GCAM"], y=df_upper[str(year)], color=colors[idx], marker="x")

        # finalize plot
        plt.ylabel(y_label)
        plt.xlabel("Region")
        plt.xticks(rotation=60, ha='right')
        plt.title(title)
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(bottom=0.5, right=.7, left=.15)
        plt.gcf().set_size_inches(12, 8)
        plt.savefig("data/data_analysis/images/" + str(RCP) + "/"  +  title + ".png", dpi=300)
        plt.show()


def plot_line_product_CI(dataframe, products, column, SSP_baseline, differentiator, title, RCP, nonBaselineScenario):
    """
    Plots a line grouped by product
    :param dataframe: the data being plotted
    :param products: the list of products in the data being plotted
    :param column: the column by which the products can be differentiated
    :param SSP_baseline: the SSPs being plotted
    :param differentiator: a secondary column to differentiate the products
    :param title: the title of the plot
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    if dataframe.empty:  # if the datafame is empty, nothing can be plotted
        print("empty dataframe")
        return
    # get plot information
    # get subplot size
    nrow, ncol = get_subplot_dimensions(["plot", "sidebar"])

    # make plots and color scheme
    fig, axs = plt.subplots(ncol, nrow, sharey='all', gridspec_kw={'wspace': 0.02, 'hspace': 0.2}, width_ratios=[12, 1])

    # find the number of model versions
    # get color scheme based on number of model versions
    versions = dataframe[differentiator].unique()
    colors, num_colors = get_colors(len(versions))
    try:
        color_counter = 0
        for i in products:
            if "SSP" in SSP_baseline:
                y = dataframe[(dataframe[column] == i) & (dataframe['SSP'] == SSP_baseline)]
                skip_years = c.GCAMConstants.skip_years
            else:
                y = dataframe[(dataframe[column] == i) & (dataframe['Version'] == SSP_baseline)]
                skip_years = c.GCAMConstants.skip_years + 1

            if not y.empty:
                # plot all versions in y
                color = colors[color_counter]

                # get line of data to plot and plot it
                y_to_plot = y.values.tolist()[0][skip_years:skip_years + len(
                    c.GCAMConstants.biochar_x)]  # only take the x values
                plot_line_on_axs(c.GCAMConstants.biochar_x, y_to_plot, str(i), color, axs, nrow, ncol, 0)

                # get min and max data across SSPs
                CI_df = dataframe[(dataframe[column] == i)]
                min_seq = [CI_df[str(i)].min() for i in c.GCAMConstants.biochar_x]
                max_seq = [CI_df[str(i)].max() for i in c.GCAMConstants.biochar_x]

                # plot the min and max data
                axs[0].fill_between(c.GCAMConstants.biochar_x, y1=min_seq, y2=max_seq, alpha=0.15, color=color)

                # add rectangles on right plot
                bar_year = str(max(c.GCAMConstants.biochar_x))
                rect = patches.Rectangle((color_counter * 0.2, CI_df[bar_year].min()),
                                         0.15, CI_df[bar_year].max() - CI_df[bar_year].min(), facecolor=color)

                # Add the patch to the Axes
                axs[1].add_patch(rect)

                color_counter = color_counter + 1

        # get units
        units = dataframe['Units'].unique()[0]
        l, h = finalize_line_subplot(axs, units, title, ncol, nrow, 0)

        axs[1].axis('off')

        finalize_line_plot(fig, h, l, axs, nrow, ncol, 2, title, RCP, nonBaselineScenario)

    except ValueError as e:
        print(e)


def plot_regional_rose(dataframe, year, SSPs, y_label, title, column, RCP, nonBaselineScenario):
    """
    Plots regional data in a categorical scatterplot
    :param dataframe: data being plotted
    :param year: evaluation column
    :param SSPs: SSPs being evaluated
    :param y_label: ylabel for graph
    :param title: title of graph
    :param column: column used to identify unique categories
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    # plot for each SSP
    for i in SSPs:
        dataframe = dataframe[dataframe['SSP'].str.contains(i)]
        for idx, item in enumerate(dataframe[column].unique()):
            df = dataframe.loc[dataframe[column] == str(item)]

            # set figure size
            fig = plt.figure()
            ax = plt.subplot(111, polar=True)
            ax.set_rlabel_position(0)
            ax.spines["polar"].set_color('#ffffff')
            # ax.set_title(str(item))
            ax.grid(color="#d5d5d5", linestyle="dashed")
            ax.text(np.radians(5), df[year].max(), y_label,
                    rotation=0, ha='center', va='center')
            plt.subplots_adjust(bottom=0.2, top=0.8)
            plt.xticks([])
            cmap = plt.colormaps.get_cmap('cool')
            # fig.suptitle(title)
            normalizer = Normalize(dataframe[year].min(), dataframe[year].max())
            im = cm.ScalarMappable(norm=normalizer, cmap=cmap)

            # Compute the angle each bar is centered on:
            heights = df[year]
            width = 2 * np.pi / len(df.index)
            indexes = list(range(1, len(df.index) + 1))
            angles = [element * width for element in indexes]

            # Draw bars
            ax.bar(
                x=angles,
                height=heights,
                width=width,
                linewidth=.5,
                edgecolor="white",
                color=im.to_rgba(heights))

            # little space between the bar and the label
            labelPadding = df[year].max() / 25

            # Add labels
            for angle, height, label in zip(angles, heights, df["GCAM"]):

                # Labels are rotated. Rotation must be specified in degrees :(
                rotation = np.rad2deg(angle)

                # Flip some labels upside down
                if np.pi / 2 <= angle < 3 * np.pi / 2:
                    alignment = "right"
                    rotation = rotation + 180
                else:
                    alignment = "left"

                # Finally add the labels
                ax.text(
                    x=angle,
                    y=height + labelPadding if height + labelPadding > labelPadding * 6 else labelPadding * 6,
                    s=label,
                    ha=alignment,
                    va='center',
                    rotation=rotation,
                    rotation_mode="anchor",
                    zorder=.2,
                    fontfamily="Arial",
                    fontstretch="extra-condensed",
                    fontsize="x-large"
                )

            plt.gcf().set_size_inches(12, 12)
            plt.savefig("data/data_analysis/images/"  + str(RCP) + "/"  +  str(item) + ".png", dpi=300)
            plt.show()


def sensitivity(dataframe, RCP, base_version, year, column, Version, nonBaselineScenario, title):
    """
    Plots a tornado plot for sensitivity analyses
    :param dataframe: dataframe to be plotted
    :param RCP: RCP pathway for baseline
    :param base_version: baseline SSP
    :param year: year for evaluation
    :param column: column containing differentiation
    :param nonBaselineScenario: the list of scenarios included in the sensitivity analysis
    :param Version: the column that contains thet version of the scenario
    :return: N/A
    """
    fig, ax = plt.subplots()

    # drop np.nan
    dataframe = dataframe[dataframe[year].notna()]

    # get base values on a per product basis
    base_vals = dataframe[dataframe[[Version]].isin([base_version]).any(axis=1)]

    # get low and high values
    # the following 2 dataframes should have the same legnth as the baes values
    low_vals = dataframe.loc[dataframe.groupby(column)[year].idxmin()]
    high_vals = dataframe.loc[dataframe.groupby(column)[year].idxmax()]
    vals = pd.merge(low_vals, high_vals, on=[column], suffixes=("_low", "_high"))
    vals = pd.merge(vals, base_vals, on=[column], suffixes=("", "_base"))
    bars = pd.DataFrame()

    # calculate bars
    bars["length"] = (vals[year + "_high"] - vals[year + "_low"])
    bars["low"] = vals[year + "_low"]
    bars["low_Version"] = vals["Version_low"]
    bars["high_Version"] = vals["Version_high"]
    bars["base"] = vals[year]
    bars["high"] = vals[year + "_high"]
    bars[column] = vals[column]
    bars["Units"] = vals["Units"]
    bars["base_unscaled"] = vals[year]
    bars = bars.dropna()
    baseline_value = 0

    # sort dataframe
    bars = bars.sort_values(by=["low", "high"], ascending=True)
    ys = range(len(bars))[::-1]  # top to bottom

    # get colormap and normalize it
    cmap = plt.colormaps.get_cmap('PiYG')
    min_low = baseline_value - bars["low"].min()
    max_high = bars["high"].max() - baseline_value
    normalizer = Normalize(-max(min_low, max_high), max(min_low, max_high))

    # Plot the bars, one by one
    for y, low, value, base, low_Version, high_Version in zip(ys, bars["low"], bars["length"], bars["base"], bars["low_Version"],
                                                      bars["high_Version"]):
        # The width of the 'low' and 'high' pieces
        low_width = base - low
        high_width = low + value - base

        # plot full colorbar so that the center of the colorbar is on the vertical line, then crop by data values
        v_offset = 0.3
        ymin = y - v_offset
        ymax = y + v_offset
        im = gradient_image(ax, direction=1,
                            extent=(
                                baseline_value - max(min_low, max_high), baseline_value + max(min_low, max_high), ymin,
                                ymax),
                            cmap=cmap,
                            cmap_range=(normalizer(-max(min_low, max_high)), normalizer(max(min_low, max_high))))
        # crop image by patch
        patch = patches.Rectangle((low, ymin), width=value, height=ymax - ymin, transform=ax.transData)
        im.set_clip_path(patch)

        # add patch for border
        edge_patch = [patches.Rectangle((low, ymin), width=value, height=ymax - ymin)]
        pc = PatchCollection(edge_patch, edgecolor="#aaaaaa", facecolors="none")
        ax.add_collection(pc)

        # Display the Version as text next to the low and high bars
        x = base - low_width - (bars["high"].max()-bars["low"].min())/ 200
        plt.text(x, y-0.12, str(low_Version), va='center', ha='right', fontsize='small')
        x = base + high_width + (bars["high"].max()-bars["low"].min())/ 200
        plt.text(x, y+0.12, str(high_Version), va='center', ha='left', fontsize='small')

        # Draw a vertical line down the middle for each segment where the baseline isn't 0
        if base != 0:
            plt.vlines(base, color='black', ymin=ymin, ymax=ymax)
            plt.text(base, y +1.4*v_offset, str(base_version), va='center', ha='center', fontsize='small')

    # if there is still a need for a singular baseline
    plt.axvline(baseline_value, color='#cccccc')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Make the y-axis display the variables
    plt.yticks(ys, bars[column])

    # Set the portion of the x- and y-axes to show
    plt.xlim(bars["low"].min() - (bars["high"].max()-bars["low"].min())/ 7, bars["high"].max() + (bars["high"].max()-bars["low"].min())/ 7)
    plt.ylim(-1, len(bars[column]))
    # plt.xlabel("change from released model in RCP " + str(RCP) + " (" + str(bars["Units"].unique()[0]) + ")")
    plt.subplots_adjust(left=.33, right=.98, bottom=.4)
    plt.savefig("data/data_analysis/images/" + str(RCP) + "/" + str(title) + ".png",
                dpi=300)
    plt.xlabel(title)
    plt.show()


def gradient_image(ax, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.
    From: https://matplotlib.org/stable/gallery/lines_bars_and_markers/gradient_bar.html

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular, *cmap*, *extent*, and *transform* may be useful.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, interpolation='bicubic', clim=(0, 1),
                   aspect='auto', **kwargs)
    return im


def plot_regional_hist_avg(prices, year, SSPs, y_label, title, column, supply, RCP, nonBaselineScenario):
    """
    Plots regional data in a stacked bar histogram
    :param prices: price data being plotted
    :param year: evaluation column
    :param SSPs: SSPs being evaluated
    :param y_label: ylabel for graph
    :param title: title of graph
    :param column: column used to identify unique categories
    :param supply: supply data being plotted for weighted averages
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: N/A
    """
    if supply == "na":
        df = prices.loc[:, [str(year), column]]
        products = df[column].unique().tolist()
        units = prices["Units"].unique()[0]
        df = df.pivot(columns=column, values=str(year))
        df = df.sort_values([column for column in df.columns])

        # get colors
        n = len(products)
        if n == 5:
            n = 1
        colors, divisions = get_colors(n)

        if (np.isnan(prices[year].max())) or (np.isnan(prices[year].min())):
            return

        # plot histogram
        bind_width = (int(prices[year].max() + .1) - int(prices[year].min() - .1)) / 40
        bins = [bind_width * i for i in
                range(int(prices[year].min() / bind_width) - 1, int(prices[year].max() / bind_width) + 2)]

        data_series = [df[column] for column in df.columns]

        plt.hist(data_series, bins, stacked=True, label=df.columns, histtype='bar',
                 color=[colors[i] for i in range(len(products))])

        # finalize plot
        plt.ylabel(y_label)
        plt.xlabel(units)
        plt.xticks(rotation=60, ha='right')
        plt.xlim(np.floor(prices[year].min() - bind_width), np.ceil(prices[year].max() + bind_width))
        plt.title(title)
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(bottom=0.4, right=.7)
        plt.savefig("data/data_analysis/images/"  + str(RCP) + "/"  +  title + ".png", dpi=300)
        plt.show()
    else:
        colors, divisions = get_colors(1)
        # plot for each SSP
        if column != "SSP":
            for k in SSPs:
                dataframe = prices[prices['SSP'].str.contains(k)]
                supply = supply[supply['SSP'].str.contains(k)]
                plot_weighted_average_hist(colors, column, dataframe, supply, title, y_label, year, RCP, nonBaselineScenario)
        else:
            plot_weighted_average_hist(colors, column, prices, supply, title, y_label, year, RCP, nonBaselineScenario)


def plot_weighted_average_hist(colors, column, dataframe, supply, title, y_label, year, RCP, nonBaselineScenario):
    """
    :param colors: plotting colors used
    :param column: column of unique histogram types
    :param dataframe: price data
    :param supply:  supply data
    :param title: graph title
    :param y_label: y-axis label
    :param year: year being evaluated
    :param RCP: the RCP pathway on which the scenarios are evaluated
    :param nonBaselineScenario: the set of scenarios in the sensitivity analysis being evaluated
    :return: histogram plot
    """
    df = dataframe.loc[:, [str(year), column]]
    products = df[column].unique().tolist()
    df = df.pivot(columns=column, values=str(year))
    # plot histogram
    bins = [50 * i for i in range(1 + int(dataframe[year].max() / 50))]
    plt.hist([df[i] for i in products], bins, stacked=True, label=df.columns, histtype='bar',
             color=[colors[i] for i in range(len(products))])
    # calculate averages
    for idx, item in enumerate(dataframe[column].unique()):
        df_price = dataframe.loc[dataframe[column] == str(item)]
        df_supply = supply.loc[supply[column] == str(item)]
        weighted_avg = pd.merge(df_price, df_supply, on=["GCAM"])
        weighted_avg[str(year)] = weighted_avg[str(year) + "_x"] * weighted_avg[str(year) + "_y"]

        print("avg:", str(item), str(weighted_avg[str(year)].sum() / weighted_avg[str(year) + "_y"].sum()) + " " +
              str(weighted_avg["Units"].unique()[0]))
    # finalize plot
    plt.ylabel("number of regions")
    plt.xlabel(y_label)
    plt.xticks(rotation=60, ha='right')
    plt.title(title)
    plt.legend(bbox_to_anchor=(1, 1))
    plt.subplots_adjust(bottom=0.4, right=.7)
    plt.savefig("data/data_analysis/images/"  + str(RCP) + "/"  +  title + ".png", dpi=300)
    plt.show()


def plot_alluvial(df, biochar_year, base_year):
    """
    plots the alluvial graph of changes to land use in ag
    :param df: dataframe being plotted
    :param biochar_year: a year in which biochar adoption is widespread
    :param base_year: a year with no biochar adoption
    :return:
    """
    df["color"] = df.apply(lambda row: data_manipulation.mgmt_to_color(row), axis=1)
    fig = px.parallel_categories(df, dimensions=[base_year, biochar_year, "Management", "Region"],
                                 labels={"Region": "Region in " + str(biochar_year),
                                         str(base_year): "Crops in "+ str(base_year),
                                         str(biochar_year): "Crops in " + str(biochar_year),
                                         "Management": "Management Type in " + str(biochar_year)},
                                 color=df["color"], width=1920)
    fig.update_layout(margin=dict(l=500, r=500, t=100, b=100),
                      font_family="Arial",
                      font_size=20
                      )
    fig.show()
