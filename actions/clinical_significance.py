import re
import io
import requests
import datetime
import time

import pandas as pd
import numpy as np
from gzip import decompress
from pathlib import Path

from st2common.runners.base_action import Action

# Palette: Brown2Blue10Steps
COLOR = {
    "Likely%20benign": "0, 169, 204",
    "Benign": "50, 227, 255",
    "Benign%2FLikely%20benign": "101, 239, 255",
    "Uncertain%20significance": "204, 253, 255",
    "Likely%20pathogenic%2C%20low%20penetrance": "242, 218, 205",
    "Pathogenic%2C%20low%20penetrance": "216, 175, 151",
    "Pathogenic%2FLikely%20pathogenic": "204, 155, 122",
    "Likely%20pathogenic": "153, 96, 53",
    "Pathogenic": "102, 47, 0",
    "Unknown": "153, 153, 153",
}


class ClinicalSignificance(Action):

    def run(
        self,
        urls: list,
        outputfolder: str,
        columns: list[str],
        attributes: list[str],
        filename: str,
    ) -> tuple[bool, str]:

        finaldf = pd.DataFrame()
        for url in urls:
            succeded, results = self.get_data(url, columns, attributes)
            if not succeded:
                return (False, results)
            if results.empty:
                continue

            df = self.filter_data(results)
            finaldf = pd.concat([finaldf, df])

            time.sleep(5)

        finaldf = self.add_comments(finaldf)
        finaldf = self.add_color(finaldf)
        filenames = self.create_annotationtrack_files(finaldf, outputfolder, filename)

        return (True, filenames)

    def get_data(
        self, url: str, columns: list, attributes: list
    ) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving data from request
        """
        try:
            response = requests.get(url)
            response.raise_for_status()
            try:
                response = decompress(response.content)
            except:
                return (False, f"File {url} is not a gzip-compressed file")
        except:
            return (False, f"Problem with request from {url}")

        # Check how many rows after comments are over, if zero then skip
        try:
            df = pd.read_csv(
                io.StringIO(response.decode()),
                sep="\t",
                header=None,
                comment="#",
                dtype=str,
            )
        except pd.errors.EmptyDataError:
            return (True, pd.DataFrame())
        except:
            return (False, f"Problem with creating dataframe of {url}")

        df.columns = columns
        df = self.get_attributes(
            df=df,
            attributes=attributes,
            last_column=columns[-1],
        )

        return (True, df)

    def get_attributes(
        self, df: pd.DataFrame, attributes: list[str], last_column: str
    ) -> pd.DataFrame:
        """
        Function for extracting key value data
        """

        for attribute in attributes:
            df[attribute] = df[last_column].apply(
                lambda x: (
                    re.findall(rf"{attribute}=([^;]*)", x)[0]
                    if f"{attribute}=" in x
                    else np.nan
                )
            )
            df[attribute] = df[attribute].astype("str")
            df.loc[df[attribute].isna(), attribute] = "-"

        df.drop(last_column, axis=1, inplace=True)

        return df

    def filter_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the dataframe
        """
        # Filter the chromosome column -> Looks like: NC_000019.10 (We want 19)
        df = df[df["chromosome"].str.contains("NC") == True].copy()
        df.loc[:, "chromosome"] = df["chromosome"].str.split(".").str[0]
        df.loc[:, "chromosome"] = df["chromosome"].str.split("_").str[1]
        df["chromosome"] = df["chromosome"].astype("int")

        # Remove any additional chromosomes
        df = df[(df["chromosome"] < 25)]

        # Change chromosome number 23 and 24 to X and Y respectively
        df["chromosome"] = df["chromosome"].astype("str")
        df.loc[df["chromosome"] == "23", "chromosome"] = "X"
        df.loc[df["chromosome"] == "24", "chromosome"] = "Y"

        df.loc[:, "Dbxref"] = df["Dbxref"].str.split(",").str[0]

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """

        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from "
            + df["source"]
            + ", see more at: https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/README.txt"
            + ";"
            + "Clinical significance: "
            + df["clinical_int"]
            + ";"
            + "SV TYPE: "
            + df["Name"]
            + ";"
            + "Phenotype: "
            + df["phenotype"]
            + ";"
            + "Phenotype ID: "
            + df["phenotype_id"]
            + ";"
            + "Consequences: "
            + df["consequence"]
            + ";"
            + "Validated: "
            + df["validated"]
            + ";"
            + "Copy number: "
            + df["copy_number"]
            + ";"
            + "Outer and/or inner range of start position: "
            + df["Start_range"]
            + ";"
            + "Outer and/or inner range of end position: "
            + df["End_range"]
            + ";"
            + "Zygosity: "
            + df["zygosity"]
            + ";"
            + "Web link to the variant: "
            + df["Dbxref"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for clinical in df["clinical_int"].unique():
            if clinical in COLOR.keys():
                df.loc[df["clinical_int"] == clinical, "color"] = COLOR[clinical]
            else:
                df.loc[df["clinical_int"] == clinical, "color"] = COLOR["Unknown"]

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> list[str]:

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        filenames = []
        for significance in df["clinical_int"].unique():

            df_copy = df[df["clinical_int"] == significance]

            df_copy = df_copy.drop(
                df_copy.columns.difference(
                    [
                        "chromosome",
                        "start",
                        "end",
                        "comments",
                        "color",
                    ]
                ),
                axis=1,
            )

            df_copy.to_csv(
                (outputfolder + f"{filename}_{significance}.tsv"),
                sep="\t",
                index=False,
            )
            filenames.append(f"{filename}_{significance}.tsv")

        return filenames
