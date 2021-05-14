# Imports
import pandas as pd
import re
import argparse

# Do not truncate long strings in pandas dataframes
pd.set_option("display.max_colwidth", 999)

# Create the parser
parser = argparse.ArgumentParser(
    description="Compares GC-content and genomes size between input genomes and NCBI databse on genus level"
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

# Parse the arguments
args = parser.parse_args()

# Open the input files
ncbi_genomes = pd.read_csv(args.ncbi_genomes, dtype={"RefSeq category": str})
genome_characteristics = pd.read_csv(args.genomes, sep="\t")
ncbi_taxonomy = pd.read_csv(args.taxonomy, sep="\t")

# Select only the genomes that are assigned at least at the genus level in the NCBI classification
mask = ncbi_taxonomy["NCBI classification"].str.split(";").str[-2].str.contains("g__$")
genus_level = ncbi_taxonomy[~mask].reset_index(drop=True)

# Remove GTDB classification
genus_level = genus_level.drop(columns=["GTDB classification"])

# Filter out irrelevant columns from the NCBI database
cols = ["#Organism Name", "Organism Groups", "Size(Mb)", "GC%", "Level"]
ncbi_genomes = ncbi_genomes[cols].copy()

# Select only complete genomes
ncbi_genomes = ncbi_genomes[ncbi_genomes["Level"] == "Complete"].copy()

# Extract genus and species in separated columns // for now only species and  genus from #OrganismName column
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

# Select relevant columns from the table of genomes under investigation
cols = [
    "Bin Id",
    "Genome size (bp)",
    "GC",
    "Completeness",
    "Contamination",
    "# scaffolds",
]
genome_characteristics = genome_characteristics[cols]

# Exclude the entries where genus.species is NaN
# This is because the Organism name is either bacterium or archaeon
ncbi_genomes = ncbi_genomes[~ncbi_genomes["genus.species"].isnull()].copy()


def remove_brackets(string):
    """Takes a string as input, removes brackets ('[]') from each end of the string"""
    return string.lstrip("[").rstrip("]")


# Remove brackets from genus, species, and genus.species columns
for col in ["genus", "species", "genus.species"]:
    ncbi_genomes[col] = ncbi_genomes[col].apply(remove_brackets).copy()

# Remove the genomes not assigned at at least genus level
genome_characteristics = genome_characteristics[
    genome_characteristics["Bin Id"].isin(genus_level["user_genome"])
]

# Merge genomes under investigation and taxnomic assignmetns
genome_characteristics = genome_characteristics.merge(
    genus_level, left_on=["Bin Id"], right_on=["user_genome"]
).drop(columns=["user_genome"])

# Add genus to user's genomes
genome_characteristics["genus"] = (
    genome_characteristics["NCBI classification"]
    .str.split(";")
    .str[-2]
    .str.split("__")
    .str[1]
)
# Extract genus and species to separate columns in genomes under investigation
genome_characteristics["genus"] = genome_characteristics["Bin Id"].str.split("_").str[0]
genome_characteristics["species"] = (
    genome_characteristics["Bin Id"].str.split("_").str[1]
)
genome_characteristics["genus.species_genomes"] = (
    genome_characteristics["genus"] + " " + genome_characteristics["species"]
)

# Set genus as index in genomes under investigation
genome_characteristics = genome_characteristics.set_index("genus")

# Group genomes by genus in genomes under investigation to compute statistics
genomes_grouped = genome_characteristics.groupby("genus")

# Compute GC content by genus in genomes under investigation and add the corresponding column
genome_characteristics["genomes_GC_genus"] = genomes_grouped["GC"].mean()

# Compute genome size by genus in genomes under investigation and add the corresponding column
genome_characteristics["genomes_size_genus"] = genomes_grouped[
    "Genome size (bp)"
].mean()

# Compute standard deviation of GC-content by genus in in genomes under investigation and add the corresponding column
genome_characteristics["genomes_GC_genus_std"] = genomes_grouped["GC"].std()

# Compute standard deviation of genome size by genus in in genomes under investigation and add the corresponding column
genome_characteristics["genomes_size_genus_std"] = genomes_grouped[
    "Genome size (bp)"
].std()

# Rename the Genome size (bp) and GC columns to distinguish the species level genome size and GC-content from genus level
genome_characteristics = genome_characteristics.rename(
    {"Genome size (bp)": "genomes_genome_size_species", "GC": "genomes_GC_species"},
    axis="columns",
)

# Set genus as index
ncbi_genomes = ncbi_genomes.set_index("genus")

# Subset from the NCBI database only the genuses from the genomes under investigation
ncbi_genomes = ncbi_genomes[
    ncbi_genomes.index.get_level_values("genus").isin(
        genome_characteristics.index.get_level_values("genus")
    )
]

# Group genomes by genus in genomes under investigation to compute statistics
ncbi_grouped = ncbi_genomes.groupby("genus")

# Compute the %GC by genus in the NCBI database
ncbi_genomes["ncbi_GC_genus"] = ncbi_grouped[
    "GC%"
].mean()  # There will be some orders/plant families but we do not care because we will exclude them when combining the data sets

# Compute genome size by genus in the NCBI database
ncbi_genomes["ncbi_genome_size_genus"] = ncbi_grouped["Size(Mb)"].mean()

# Compute standard deviation of GC-content by genus in NCBI database and add the corresponding column
ncbi_genomes["ncbi_genome_GC_genus_std"] = ncbi_grouped["GC%"].std()

# Compute standard deviation of genome size by genus in NCBI database and add the corresponding column
ncbi_genomes["ncbi_genome_size_genus_std"] = ncbi_grouped["Size(Mb)"].std()

# Rename the Size(Mb) and GC(%) columns to distinguish the species level genome size and GC-content from genus level
ncbi_genomes = ncbi_genomes.rename(
    {"Size(Mb)": "ncbi_genome_size_species", "GC%": "ncbi_GC_species"}, axis="columns"
)

# Merge NCBI database and genomes under investigation
merged = pd.merge(
    genome_characteristics, ncbi_genomes, left_index=True, right_index=True
)

# Drop species duplicates
merged = merged.drop_duplicates(subset=["Bin Id"], keep="first")


def report():
    """Saves the report on GC-content and genome size in a csv file 'report.csv'"""
    # Columns to include in the report
    report_cols = [
        "Completeness",
        "Contamination",
        "# scaffolds",
    ]

    df = merged[report_cols].copy()

    df["GC_diff"] = merged["genomes_GC_genus"] - merged["ncbi_GC_genus"]
    df["genome_size_diff"] = (
        merged["genomes_size_genus"] - merged["ncbi_genome_size_genus"]
    ).astype("int64")
    df["GC_std"] = (df["GC_diff"] / merged["ncbi_genome_GC_genus_std"]).abs()
    df["genome_size_std"] = (
        df["genome_size_diff"] / merged["ncbi_genome_size_genus_std"]
    ).abs()

    return df.to_csv("report.csv")


if __name__ == "__main__":
    report()
