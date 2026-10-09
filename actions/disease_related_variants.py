import io
import requests
import datetime

import pandas as pd
import numpy as np
from pathlib import Path

from st2common.runners.base_action import Action

CHROMS = {
    "1": 248_956_422,
    "2": 242_193_529,
    "3": 198_295_559,
    "4": 190_214_555,
    "5": 181_538_259,
    "6": 170_805_979,
    "7": 159_345_973,
    "8": 145_138_636,
    "9": 138_394_717,
    "10": 133_797_422,
    "11": 135_086_622,
    "12": 133_275_309,
    "13": 114_364_328,
    "14": 107_043_718,
    "15": 101_991_189,
    "16": 90_338_345,
    "17": 83_257_441,
    "18": 80_373_285,
    "19": 58_617_616,
    "20": 64_444_167,
    "21": 46_709_983,
    "22": 50_818_468,
    "X": 156_040_895,
    "Y": 57_227_415,
}

COLOR = {
    "benign": "50, 227, 255",
    "variants of uncertain significance": "204, 253, 255",
    "na; variants of uncertain significance; na": "204, 253, 255",
    "na; pathogenic; variants of uncertain significance": "204, 253, 255",
    "pathogenic/likely pathogenic": "204, 155, 122",
    "likely pathogenic": "153, 96, 53",
    "pathogenic": "102, 47, 0",
    "Unknown": "153, 153, 153",
}


class DiseaseRelatedSVs(Action):

    def run(self, url: str, outputfolder: str, filename: str) -> tuple[bool, str]:

        succeded, results = self.get_data(url)
        if not succeded:
            return (False, results)

        df = self.get_genomic_position(results)
        df = self.add_comments(df)
        df = self.add_color(df)
        filename = self.create_annotationtrack_files(df, outputfolder, filename)

        return (True, filename)

    def get_data(self, url: str) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving data from request
        """
        # https://www.geeksforgeeks.org/python/exception-handling-of-python-requests-module/
        try:
            response = requests.get(url)
            response.raise_for_status()
        except requests.exceptions.HTTPError as error:
            (
                False,
                f"HTTP problem with request from {url}, errormessage: {error.args[0]}",
            )
        except requests.exceptions.ReadTimeout as error:
            (False, f"Timeout for {url}, errormessage: {error.args[0]}")
        except requests.exceptions.ConnectionError:
            (False, f"Problem with internet connection during request of {url}")
        except:
            return (False, f"Problem with request from {url}")

        df = pd.json_normalize(response.json())

        return (True, df)

    def get_genomic_position(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for extracting the genomic positions for the genes
        """
        df = df[df["Genome_Coordinate_Hg38"].notna()].copy()

        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "chromosome", ";", ":", "chr", 0, 0, 1
        )
        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "start", ";", ":", "-", 0, 1, 0
        )
        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "end", ";", ":", "-", 0, 1, 1
        )
        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "chromosomeEND", ";", ":", "chr", 1, 0, 1
        )
        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "startEND", ";", ":", "-", 1, 1, 0
        )
        df = self.split_dataframe(
            df, "Genome_Coordinate_Hg38", "endEND", ";", ":", "-", 1, 1, 1
        )

        # Add an start/end positions to those SVs that don't have one
        df.loc[((df["endEND"].isna()) & (df["startEND"].notna())), "endEND"] = df[
            "startEND"
        ]
        df.loc[((df["end"].isna()) & (df["start"].notna())), "end"] = df["start"]

        # df.loc[(df["end"] == "qter"), "end"] = CHROMS[df["chromosome"]]
        df.loc[(df["end"] == "qter"), "end"] = df.loc[
            df["end"] == "qter", "chromosome"
        ].map(CHROMS)

        # If there is no end chromosome, set the end chromosome value to the start chromosome
        df.loc[:, "chromosomeEND"] = (
            df["chromosomeEND"].str.strip().replace("NA", np.nan)
        )
        df.loc[(df["chromosomeEND"].isna()), "chromosomeEND"] = df["chromosome"]

        # Set True or False if the SV starts and stops onn the same chromosome
        df["interchromosomal"] = np.where(
            (df["chromosome"] == df["chromosomeEND"]), "True", "False"
        )

        return df

    def split_dataframe(
        self,
        df: pd.DataFrame,
        inputcolumn: str,
        outputcolumn: str,
        firstsplit: str,
        secondsplit: str,
        thirdsplit: str,
        firstvalue: int,
        secondvalue: int,
        thirdvalue: int,
    ) -> pd.DataFrame:

        df.loc[:, outputcolumn] = (
            df[inputcolumn]
            .str.split(firstsplit)
            .str[firstvalue]
            .str.split(secondsplit)
            .str[secondvalue]
            .str.split(thirdsplit)
            .str[thirdvalue]
        )

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from SV4GD, see more at: https://bio-computing.hrbmu.edu.cn/SV4GD/download"
            + ";"
            + "Variant: "
            + df["Type"]
            + ": "
            + df["Variant_Type"]
            + ";"
            + "Disease Category: "
            + df["Disease_Category"]
            + ";"
            + "Pathogenicity: "
            + df["Pathogenicity"]
            + ";"
            + "ACMG Score (Standard): "
            + df["ACMG Score (Standard)"]
            + ";"
            + "Main Phenotypes: "
            + df["Main_Phenotypes"]
            + ";"
            + "Disease Diagnosis: "
            + df["Disease_Diagnosis"]
            + ";"
            + "Description of pathogenic mechanism [Most_Susceptible_gene]: "
            + df["Description_of_pathogenic_mechanism [Most_Susceptible_gene]"]
            + ";"
            + "Country/Region: "
            + df["Country_Region"]
            + ";"
            + "PMID: "
            + df["PMID"]
            + ";"
            + "Article info: "
            + ": "
            + df["Journal"]
            + ": "
            + df["Article_Title"]
            + ": "
            + df["Year"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for significance in df["Pathogenicity"].unique():
            if significance.lower() in COLOR.keys():
                df.loc[df["Pathogenicity"] == significance, "color"] = COLOR[
                    significance.lower()
                ]
            else:
                df.loc[df["Pathogenicity"] == significance, "color"] = COLOR["Unknown"]

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> str:

        df_interchromosomal = df[df["interchromosomal"] == "True"]
        df_non_interchromosomal_start = df[df["interchromosomal"] == "False"]
        df_non_interchromosomal_end = df[df["interchromosomal"] == "False"]

        # Overwrite values
        df_non_interchromosomal_end.loc[:, "start"] = df_non_interchromosomal_end[
            "startEND"
        ]
        df_non_interchromosomal_end.loc[:, "end"] = df_non_interchromosomal_end[
            "endEND"
        ]
        df_non_interchromosomal_end.loc[:, "chromosome"] = df_non_interchromosomal_end[
            "chromosomeEND"
        ]

        df = pd.concat(
            [
                df_interchromosomal,
                df_non_interchromosomal_start,
                df_non_interchromosomal_end,
            ]
        )

        df = df.drop(
            df.columns.difference(
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

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df.to_csv(
            (outputfolder + filename),
            sep="\t",
            index=False,
        )

        return [filename]
