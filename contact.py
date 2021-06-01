# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
from pathlib import Path
from os.path import splitext

# Create the parser
parser = argparse.ArgumentParser(
    description="Compares GC-content and genomes size between input genomes and reference genomes from NCBI database on genus level. Outputs results in tabular and graphic format."
)

# Add the arguments
parser.add_argument("--genomes", help="Input genomes to process.", required=True)
parser.add_argument(
    "--ncbi_genomes", help="Input NCBI reference genomes.", required=True
)
parser.add_argument("--taxonomy", help="Input taxonomy file.", required=True)
parser.add_argument(
    "--n_entries", help="Threshold for number of entries in NCBI, default 5.", default=5
)
parser.add_argument(
    "--completeness",
    help="Pass 'all' if you want to use all NCBI genomes and 'complete' if you want to use only complete genomes, default 'complete'.",
    default="complete",
    choices=["complete", "all"],
)

# Parse the arguments
args = parser.parse_args()

# Remove first 6 rows of technical information from checkM output file (input genomes)
# This will create a new file with '_table' added to its name
try:
    with open(args.genomes, "r") as checkM_old:
        data = checkM_old.read().splitlines(True)

    with open(splitext(args.genomes)[0] + "_table.txt", "w") as checkM:
        checkM.writelines(data[6:])

    input_genomes = pd.read_csv(splitext(args.genomes)[0] + "_table.txt", sep="\t")

except FileNotFoundError:
    print(f"Input genomes file {args.genomes} not found. Check if the path is correct.")
    exit()

# Open the input files
try:
    ncbi_genomes = pd.read_csv(args.ncbi_genomes, dtype={"RefSeq category": str})
except FileNotFoundError:
    print(
        f"Input NCBI genomes file {args.ncbi_genomes} not found. Check if the path is correct."
    )
    exit()

try:
    ncbi_taxonomy = pd.read_csv(args.taxonomy, sep="\t")
except FileNotFoundError:
    print(
        f"Input taxonomy file {args.taxonomy} not found. Check if the path is correct."
    )
    exit()

# Select only the genomes that are assigned at least at the genus level in the NCBI classification
mask = ncbi_taxonomy["NCBI classification"].str.split(";").str[-2].str.contains("g__$")
genus_level = ncbi_taxonomy[~mask].reset_index(drop=True)

# Remove GTDB classification
genus_level = genus_level.drop(columns=["GTDB classification"])

# Filter out irrelevant columns from the NCBI database
cols = ["#Organism Name", "Organism Groups", "Size(Mb)", "GC%", "Level"]
ncbi_genomes = ncbi_genomes[cols].copy()

# Select only complete genomes if specified in argument 'completeness'
if args.completeness == "complete":
    ncbi_genomes = ncbi_genomes[ncbi_genomes["Level"] == "Complete"].copy()

# Extract genus and species in separated columns // for now only species and genus from #OrganismName column
ncbi_genomes["Organism Groups"].str.split(";").str[0]
ncbi_genomes["genus"] = ncbi_genomes["#Organism Name"].str.split(" ").str[0]
ncbi_genomes["species"] = ncbi_genomes["#Organism Name"].str.split(" ").str[1]
ncbi_genomes["genus.species"] = ncbi_genomes["genus"] + " " + ncbi_genomes["species"]

# Remove uncultured bacteria
ncbi_genomes = ncbi_genomes[~(ncbi_genomes["genus"] == "uncultured")]

# Transform genome size in millions
ncbi_genomes["Size(Mb)"] = ncbi_genomes["Size(Mb)"] * 10 ** 6

# Exclude genomes that have number of entries less than the input threshold
ncbi_genomes = ncbi_genomes.groupby("genus").filter(
    lambda x: len(x) >= int(args.n_entries)
)

# Select relevant columns from the table of input genomes
cols = [
    "Bin Id",
    "Genome size (bp)",
    "GC",
    "Completeness",
    "Contamination",
    "# scaffolds",
]
input_genomes = input_genomes[cols]

# Exclude the entries where genus.species is NaN
# This is because the Organism name is either bacterium or archaeon
ncbi_genomes = ncbi_genomes[~ncbi_genomes["genus.species"].isnull()].copy()


def remove_brackets(string):
    """
    Takes a string as input, removes square brackets ('[]') from each side of the string.
    """
    return string.lstrip("[").rstrip("]")


# Remove brackets from genus, species, and genus.species columns
for col in ["genus", "species", "genus.species"]:
    ncbi_genomes[col] = ncbi_genomes[col].apply(remove_brackets).copy()

# Remove the genomes not assigned at at least genus level
input_genomes = input_genomes[input_genomes["Bin Id"].isin(genus_level["user_genome"])]

# Merge input genomes and taxonomic assignments
input_genomes = input_genomes.merge(
    genus_level, left_on=["Bin Id"], right_on=["user_genome"]
).drop(columns=["user_genome"])

# Add genus column to input genomes
input_genomes["genus"] = (
    input_genomes["NCBI classification"].str.split(";").str[-2].str.split("__").str[1]
)

# Add species column to input genomes
input_genomes["species"] = (
    input_genomes["NCBI classification"].str.split(";").str[-1].str.split("__").str[1]
)

# Species is sp if no species is provided
input_genomes["species"].apply(lambda x: "sp" if x == "" else x)

# Add genus.species column to input genomes
input_genomes["genus.species_genomes"] = (
    input_genomes["genus"] + " " + input_genomes["species"]
)

# Set genus as index in input genomes
input_genomes = input_genomes.set_index("genus")

# Rename the Genome size (bp) and GC columns to distinguish the species level genome size and GC-content from genus level
input_genomes = input_genomes.rename(
    {"Genome size (bp)": "genomes_genome_size_species", "GC": "genomes_GC_species"},
    axis="columns",
)

# Set genus as index
ncbi_genomes = ncbi_genomes.set_index("genus")

# Subset from the NCBI database only the genera present in input genomes
ncbi_genomes = ncbi_genomes[
    ncbi_genomes.index.get_level_values("genus").isin(
        input_genomes.index.get_level_values("genus")
    )
]

# Group genomes by genus in input genomes to compute downstream statistics
ncbi_grouped = ncbi_genomes.groupby("genus")

# Compute the %GC by genus in the NCBI database and add the corresponding column
ncbi_genomes["ncbi_GC_genus"] = ncbi_grouped[
    "GC%"
].mean()  # There will be some plant families but we do not care because we will exclude them when combining the data sets

# Compute mean genome size by genus in the NCBI database and add the corresponding column
ncbi_genomes["ncbi_genome_size_genus"] = ncbi_grouped["Size(Mb)"].mean()

# Compute standard deviation of GC-content by genus in NCBI database and add the corresponding column
ncbi_genomes["ncbi_genome_GC_genus_std"] = ncbi_grouped["GC%"].std()

# Compute standard deviation of genome size by genus in NCBI database and add the corresponding column
ncbi_genomes["ncbi_genome_size_genus_std"] = ncbi_grouped["Size(Mb)"].std()

# Merge NCBI database and input genomes
merged = pd.merge(input_genomes, ncbi_genomes, left_index=True, right_index=True)

# Drop species duplicates
merged = merged.drop_duplicates(subset=["Bin Id"], keep="first")


def report():
    """
    Saves the report on GC-content and genome size in a tsv file 'report.tsv'.
    """
    # Columns to include in the report
    report_cols = [
        "Bin Id",
        "Completeness",
        "Contamination",
        "# scaffolds",
    ]

    df = merged[report_cols].copy()

    # Compute statistics
    df["GC_diff"] = (merged["genomes_GC_species"] - merged["ncbi_GC_genus"]).round(3)
    df["genome_size_diff"] = (
        merged["genomes_genome_size_species"] - merged["ncbi_genome_size_genus"]
    ).astype("int64")
    df["GC_std"] = (df["GC_diff"] / merged["ncbi_genome_GC_genus_std"]).abs().round(3)
    df["genome_size_std"] = (
        (df["genome_size_diff"] / merged["ncbi_genome_size_genus_std"]).abs().round(3)
    )

    # If we divide by 0 in pandas dataframes, we get +-inf that we have to convert in NaN values
    df.loc[~np.isfinite(df["GC_std"]), "GC_std"] = np.nan
    df.loc[~np.isfinite(df["genome_size_std"]), "genome_size_std"] == np.nan

    return df


def plot_style(ax):
    """
    Defines boxplot style.
    """
    # Remove ticks and set label size
    ax.tick_params(axis="both", labelsize=14, bottom=False, left=False)

    # Set x axis label size
    ax.set_xlabel(xlabel=plt.gca().get_xlabel(), size=18)

    ax.set_ylabel(ylabel=plt.gca().get_ylabel(), size=18)

    # Remove spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def draw_boxplot(col, genus, bin_id):
    """
    Creates a folder and saves boxplots for the specified measure (GC-content or genome size) of the NCBI database with a horizontal line
    that shows the position of input genome compared to the distribution of a measure of its genome in the NCBI database.

            Parameters:
                    col (str): a column, either 'GC%' for GC-content or 'Size(Mb)' for genomes size
                    genus (str): genus from NCBI database for which the distribution should be plotted
                    bin_id (str): genome id from input genomes

            Returns:
                    Creates folders for boxplots if do not exist and returns Matplotlib Axes
    """

    # Create folders for boxplots if do not exist
    if args.completeness == "complete":
        if col == "GC%":
            Path("./boxplots_complete/GC-content").mkdir(parents=True, exist_ok=True)
        else:
            Path("./boxplots_complete/genome-size").mkdir(parents=True, exist_ok=True)
    else:
        if col == "GC%":
            Path("./boxplots_all/GC-content").mkdir(parents=True, exist_ok=True)
        else:
            Path("./boxplots_all/genome-size").mkdir(parents=True, exist_ok=True)

    # Boxplots for each genus from NCBI database
    # If the genus was filtered out we do not need its boxplot
    ax = ncbi_genomes[ncbi_genomes.index.get_level_values("genus") == genus].boxplot(
        column=col, grid=False, figsize=(8, 5)
    )
    ax.set_xticklabels("")

    plot_style(ax)
    ax.set_title(bin_id, fontdict={"size": 22, "weight": "bold", "alpha": 0.75})

    if col == "GC%":
        ax.set_ylabel("GC-content (%)")

        # Add horizontal lines from input genomes for each genus to visually show the difference between
        # input and database
        ax.axhline(
            input_genomes[input_genomes["Bin Id"] == bin_id]["genomes_GC_species"]
            .astype(int)
            .values,
            color="red",
        )

    else:
        ax.axhline(
            input_genomes[input_genomes["Bin Id"] == bin_id][
                "genomes_genome_size_species"
            ]
            .astype(int)
            .values,
            color="red",
        )
    return ax


def main():
    # Save the report as a tsv file
    report().to_csv("report.tsv", sep="\t")

    # Save boxplots of GC% for complete genomes and output in console excluded input ids
    if args.completeness == "complete":
        for bin_id in input_genomes["Bin Id"]:
            genus = input_genomes[input_genomes["Bin Id"] == bin_id].index.values[0]
            if genus in ncbi_genomes.index.get_level_values("genus").unique():
                draw_boxplot("GC%", genus, bin_id)
                plt.savefig(f"./boxplots_complete/GC-content/{bin_id}.jpg")
                plt.close()
            else:
                print(
                    f"Excluding {bin_id}, no genus found for this id in filtered NCBI database."
                )

        # Save boxplots of genomes size for complete genomes
        for bin_id in input_genomes["Bin Id"]:
            genus = input_genomes[input_genomes["Bin Id"] == bin_id].index.values[0]
            if genus in ncbi_genomes.index.get_level_values("genus").unique():
                draw_boxplot("Size(Mb)", genus, bin_id)
                plt.savefig(f"./boxplots_complete/genome-size/{bin_id}.jpg")
                plt.close()

    else:
        # Save boxplots of GC% for all genomes and output in console excluded input ids
        for bin_id in input_genomes["Bin Id"]:
            genus = input_genomes[input_genomes["Bin Id"] == bin_id].index.values[0]
            if genus in ncbi_genomes.index.get_level_values("genus").unique():
                draw_boxplot("GC%", genus, bin_id)
                plt.savefig(f"./boxplots_all/GC-content/{bin_id}.jpg")
                plt.close()
            else:
                print(
                    f"Excluding {bin_id}, no genus found for this id in filtered NCBI database."
                )
        # Save boxplots of genomes size for all genomes
        for bin_id in input_genomes["Bin Id"]:
            genus = input_genomes[input_genomes["Bin Id"] == bin_id].index.values[0]
            if genus in ncbi_genomes.index.get_level_values("genus").unique():
                draw_boxplot("Size(Mb)", genus, bin_id)
                plt.savefig(f"./boxplots_all/genome-size/{bin_id}.jpg")
                plt.close()


if __name__ == "__main__":
    main()
