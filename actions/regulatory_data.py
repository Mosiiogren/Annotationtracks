import re
import io
import datetime
import requests

import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path
from st2common.runners.base_action import Action

# colorBlindness::Blue2Gray8Steps
COLOR = {
    "promoter": "121, 130, 52",
    "enhancer": "163, 173, 98",
    "tf_binding_site": "240, 198, 195",
    "ctcf_binding_site": "208, 211, 162",
    "open_chromatin_region": "223, 145, 163",
    "emar": "212, 103, 128",
    "Unknown": "153, 153, 153",
}


class RegulatoryData(Action):

    def run(
        self,
        regulatoryfileurl: list[str],
        outputfileregulatory: str,
        outputfoldertracks: str,
        columns: list[str],
        attributes: list[str],
        filename: str,
    ) -> tuple[bool, str]:

        df = pd.DataFrame()
        for url in regulatoryfileurl:
            succeded, results = self.get_data(url, columns, attributes)
            if not succeded:
                return (False, results)

            df = pd.concat([df, results])

        df = self.add_comments(df)
        df = self.add_color(df)
        self.save_to_file(df, outputfileregulatory)
        filename = self.create_annotationtrack_files(
            df,
            outputfoldertracks,
            filename,
        )

        return (True, filename)

    def get_data(
        self, url: str, columns: list[str], attributes: list[str]
    ) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving data from request
        """
        try:
            response = requests.get(url)
            if response.status_code == 429:
                return (False, f"Too many requests! {response.status_code}")

            response.raise_for_status()
            try:
                response = decompress(response.content)
            except:
                return (False, f"File {url} is not a gzip-compressed file")
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

        df = pd.read_csv(
            io.StringIO(response.decode()),
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
        )

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
                    else "-"
                )
            )
            df[attribute] = df[attribute].astype("str")
            df.loc[df[attribute].isna(), attribute] = "-"

        df.drop(last_column, axis=1, inplace=True)

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from ENSEMBL, see more at: https://ftp.ensembl.org/pub/release-116/regulation/README.md"
            + ";"
            + "Regulatory element: "
            + df["Type"]
            + ";"
            + "ID: "
            + df["ID"]
            + ";"
            + "Extended start: "
            + df["extended_start"]
            + ";"
            + "Extended_end: "
            + df["extended_end"]
            + ";"
            + "Gene ID: "
            + df["gene_id"]
            + ";"
            + "Gene name: "
            + df["gene_name"]
            + ";"
            + "Gene biotype: "
            + df["gene_biotype"]
            + ";"
            + "Source: "
            + df["source"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for element in df["Type"].unique():
            if element.lower() in COLOR.keys():
                df.loc[df["Type"] == element, "color"] = COLOR[element.lower()]
            else:
                df.loc[df["Type"] == element, "color"] = COLOR["Unknown"]

        return df

    def save_to_file(self, df: pd.DataFrame, outputfile):
        """
        Function for saving data to json file
        """
        df = df[(df["Type"] != "open_chromatin_region") & (df["Type"] != "EMAR")]
        df.to_json(outputfile, orient="records")

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
                    "color",
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
