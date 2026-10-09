import io
import requests
import datetime
import time

import pandas as pd
from gzip import decompress

from st2common.runners.base_action import Action

COLOR = {
    "DEL": "230, 159, 0",
    "DUP": "86, 180, 233",
    "INS": "0, 158, 115",
    "Unknown": "153, 153, 153",
}


class PublicStructuralVariantsData(Action):

    def run(
        self,
        urls: list[str],
        outputfolder: str,
        columns: list[str],
        filename: str,
    ) -> tuple[bool, str]:

        finaldf = pd.DataFrame()
        for url in urls:
            succeded, results = self.get_data(url, columns)
            if not succeded:
                return (False, results)
            if results.empty:
                continue

            df = self.filter_data(results)
            df = self.add_clinical_significance(df, url)

            finaldf = pd.concat([finaldf, df])

            time.sleep(5)

        df = self.add_color(finaldf)
        df = self.add_comments(df)
        filenames = self.create_annotationtrack_files(df, outputfolder, filename)

        return (True, filenames)

    def get_data(self, url: str, columns: list[str]) -> tuple[bool, pd.DataFrame | str]:
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

        return (True, df)

    def filter_data(self, df: pd.DataFrame) -> pd.DataFrame:

        df.loc[:, "chromosome"] = df["chromosome"].str.split("chr").str[1]
        df["Name"] = df["info"].str.split("_").str[-1]

        return df

    def add_clinical_significance(self, df: pd.DataFrame, url: str) -> pd.DataFrame:

        df["clincal_significance"] = (
            url.rsplit("/")[-1].rsplit(".bed")[0].rsplit(".")[-1]
        )

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """

        df = df.astype(str)

        df["comments"] = (
            "Data retrievede from dbVar, see more at: https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/annotation/GRCh38/README"
            + ";"
            + "SV TYPE: "
            + df["Name"]
            + ";"
            + "Clinical significance: "
            + df["clincal_significance"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
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
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> list[str]:

        filenames = []
        for SVtype in df["Name"].unique():
            for significance in df["clincal_significance"].unique():
                df_copy = df[
                    (df["Name"] == SVtype)
                    & (df["clincal_significance"] == significance)
                ]

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
                    (outputfolder + f"{filename}_{SVtype}_{significance}.tsv"),
                    sep="\t",
                    index=False,
                )
                filenames.append(f"{filename}_{SVtype}_{significance}.tsv")

        return filenames
