import re
import io
import requests
import datetime

import pandas as pd
import numpy as np
from gzip import decompress
from pathlib import Path

from st2common.runners.base_action import Action


class AcmgGenes(Action):

    def run(
        self,
        url: str,
        outputfolder: str,
        columns: list[str],
        attributes: list[str],
        filename: str,
    ) -> tuple[bool, str]:

        succeded, results = self.get_data(url, columns, attributes)
        if not succeded:
            return (False, results)

        finaldf = self.filter_data(results)
        finaldf = self.add_comments(finaldf)
        filename = self.create_annotationtrack_files(finaldf, outputfolder, filename)

        return (True, filename)

    def get_data(
        self, url: str, columns: list[str], attributes: list[str]
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

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from dbVar, see more at: https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/annotation/GRCh38/README"
            + ";"
            + "Gene Name: "
            + df["gene"]
            + ";"
            + "Gene synonym: "
            + df["gene_synonym"]
            + ";"
            + "Gene type: "
            + df["gene_biotype"]
            + ";"
            + "Gene description: "
            + df["description"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> list[str]:

        df = df.drop(
            df.columns.difference(
                [
                    "chromosome",
                    "start",
                    "end",
                    "comments",
                ]
            ),
            axis=1,
        )

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df.to_csv(
            (outputfolder + filename),
            sep="\t",
            index=False,
        )

        return [filename]
