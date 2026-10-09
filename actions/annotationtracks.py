import datetime
import pandas as pd
import numpy as np

from st2common.runners.base_action import Action

# colorblindr::Okabelto
COLOR = {
    "DEL": "230, 159, 0",
    "DUP": "86, 180, 233",
    "INS": "0, 158, 115",
    "DUP:TANDEM": "240, 228, 66",
    "INV": "0, 114, 178",
    "BND": "213, 94, 0",
    "DUP:INV": "204, 121, 167",
    "Unknown": "153, 153, 153",
}


class Annotationtracks(Action):
    """
    Class for creating annotationtracks from retrieved clustered data
    """

    def run(self, clusteringdata, outputfolder) -> tuple[bool, list[str]]:

        df = pd.read_json(clusteringdata, orient="records")
        df["start"] = df["start"].astype("float").astype("int")
        df["end"] = df["end"].astype("float").astype("int")

        df = self.add_interchromosomal(df)
        df = self.add_comments(df)
        df = self.add_color(df)

        filenames = self.create_annotationtrack_files(df, outputfolder)

        return (True, filenames)

    def add_interchromosomal(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding interchromosomal information -> important for BND where start and stop chromosomes are not the same
        """

        df_intrachromosomal = df[df["intrachromosomal"] == "True"]
        df_non_intrachromosomal = df[df["intrachromosomal"] == "False"]

        # BND cannot be placed in GENS, however one can add a cluster that starts at one chromosome
        # And another cluster that ends at another chromosome to be able to visualize them in Gens
        df_non_intrachromosomal_copy = df_non_intrachromosomal.copy()
        df_non_intrachromosomal_copy = df_non_intrachromosomal_copy.drop("end", axis=1)
        df_non_intrachromosomal_copy["end"] = df_non_intrachromosomal_copy["start"] + 1
        df_non_intrachromosomal_copy["chromosomeEND"] = df_non_intrachromosomal_copy[
            "chromosome"
        ]

        df_non_intrachromosomal.loc[:, "start"] = df_non_intrachromosomal["end"] - 1
        df_non_intrachromosomal["chromosome"] = df_non_intrachromosomal["chromosomeEND"]

        df = pd.concat(
            [df_intrachromosomal, df_non_intrachromosomal, df_non_intrachromosomal_copy]
        )

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """

        df = df.astype(str)
        df["comments"] = np.where(
            (df["intrachromosomal"] == "True"),
            "Data retrieved from Local Database at Klinisk Genetik"
            + ";"
            + "SVs included in the cluster: "
            + df["score"].astype(str)
            + ";"
            + "SV type: "
            + df["Name"].astype(str)
            + ";"
            + "Overlapping genes: "
            + df["genes"].astype(str)
            + ";"
            + "Overlapping exons: "
            + df["exons"].astype(str)
            + ";"
            + "Overlapping introns: "
            + df["introns"].astype(str)
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c"),
            "Data retrieved from Local Database at Klinisk Genetik"
            + ";"
            + "SVs included in the cluster: "
            + df["score"].astype(str)
            + ";"
            + "SV type: "
            + df["Name"].astype(str)
            + ";"
            + "Overlapping genes: "
            + df["genes"].astype(str)
            + ";"
            + "Overlapping exons: "
            + df["exons"].astype(str)
            + ";"
            + "Overlapping introns: "
            + df["introns"].astype(str)
            + ";"
            + "Start Chromosome: "
            + df["chromosome"].astype(str)
            + ":"
            + df["start"].astype(str)
            + ";"
            + "End chromosome: "
            + df["chromosomeEND"]
            + ":"
            + df["end"].astype(str)
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c"),
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for SVtype in df["Name"].unique():
            if SVtype.upper() in COLOR.keys():
                df.loc[df["Name"] == SVtype, "color"] = COLOR[SVtype.upper()]
            else:
                df.loc[df["Name"] == SVtype, "color"] = COLOR["Unknown"]

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str
    ) -> list[str]:
        """
        Function for creating tab separate annotation files to be uploaded into Gens
        One file for each SV type is created
        """
        # Remove unwanted chromsomes which cannot be displaced in Gens
        df = df[df["chromosome"].str.contains("Un|EBV|random|M") == False].copy()

        df["start"] = df["start"].astype("float").astype("int")
        df["end"] = df["end"].astype("float").astype("int")

        filenames = []
        for category in df["category"].unique():
            for SVtype in df["Name"].unique():
                df_copy = df[(df["Name"] == SVtype) & (df["category"] == category)]

                if df_copy.empty:
                    continue

                df_copy = df_copy.drop(
                    df_copy.columns.difference(
                        [
                            "chromosome",
                            "start",
                            "end",
                            "color",
                            "comments",
                        ]
                    ),
                    axis=1,
                )
                df_copy.to_csv(
                    (outputfolder + f"Local_{category}_{SVtype}.tsv"),
                    sep="\t",
                    index=False,
                )
                filenames.append(f"Local_{category}_{SVtype}.tsv")

        return filenames
