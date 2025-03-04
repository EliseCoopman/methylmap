import sys
import logging
import pandas as pd
import plotly.figure_factory as ff


def make_dendro(methfreqtable):
    methfreqtable = methfreqtable.apply(pd.to_numeric, errors="coerce")

    #remove columns (individuals) with 40% or more missing data
    threshold = 0.4
    missing_data_percentage = methfreqtable.isna().mean()
    columns_to_drop = missing_data_percentage[
        missing_data_percentage >= threshold
    ].index.tolist()
    methfreqtable.drop(columns=columns_to_drop, inplace=True)

    number_of_nan_values = methfreqtable.isna().sum().sum()
    if number_of_nan_values != 0:
        logging.warning(
            f"\n\n{number_of_nan_values} NaN values found in data. NaN values will be estimated using numpy interpolate to perform hierarchical clustering.\n\n"
        )
        sys.stderr.write(
            f"\n\n{number_of_nan_values} NaN values found in data. NaN values will be estimated using numpy interpolate to perform hierarchical clustering.\n\n"
        )
        methfreqtable.interpolate(method="linear", axis=0, inplace=True)
    number_of_nan_values_interpolate = methfreqtable.isna().sum().sum()
    if number_of_nan_values_interpolate != 0:
        logging.warning(
            f"\n\n{number_of_nan_values_interpolate} NaN values found in data after using numpy interpolate for estimation of these values. Individuals with minimal 1 NaN value will be deleted to perform hierarchical clustering.\n\n"
        )
        sys.stderr.write(
            f"\n\n{number_of_nan_values_interpolate} NaN values found in data after using numpy interpolate for estimation of these values. Individuals with minimal 1 NaN value will be deleted to perform hierarchical clustering.\n\n"
        )
    methfreqtable.dropna(axis=1, inplace=True)
    methfreqtable_transposed = methfreqtable.transpose()
    samples = methfreqtable_transposed.index.tolist()
    methfreqtable_transposed.reset_index(drop=True, inplace=True)

    den = ff.create_dendrogram(methfreqtable_transposed, labels=samples)

    list_sorted_samples = den.layout.xaxis.ticktext.tolist()

    methfreqtable = methfreqtable.reindex(columns=list_sorted_samples)

    return methfreqtable, den, list_sorted_samples
