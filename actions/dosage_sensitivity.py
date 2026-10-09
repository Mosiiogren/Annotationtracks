import io
import requests
import datetime

import pandas as pd
from pathlib import Path

from st2common.runners.base_action import Action

COLOR = {
    "Triplosensitivity": "80, 255, 80",
    "Haploinsufficiency": "241, 0, 241",
    "Triplosensitivity&Haploinsufficiency": "80, 1, 80",
    "Dosage sensitivity unlikely": "0, 80, 0",
    "Unknown": "100, 100, 150",
}


class DosageSensativeData(Action):

    def run(self, url: str, outputfolder: str, filename: str) -> tuple[bool, str]:

        succeded, results = self.get_data(url)
        if not succeded:
            return (False, results)

        df = self.get_genomic_position(results)
        df = self.combine_PMID(df)
        df = self.add_color(df)
        df = self.add_comments(df)
        filename = self.create_annotationtrack_files(df, outputfolder, filename)

        return (True, filename)

    def get_data(self, url: str) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving data from request
        """
        try:
            response = requests.get(url)
            response.raise_for_status()

        except:
            return (False, f"Problem with request from {url}")

        df = pd.read_csv(
            io.StringIO(response.content.decode()),
            sep="\t",
            header=0,
            skiprows=5,
            dtype=str,
        )

        return (True, df)

    def get_genomic_position(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for extracting the genomic positions for the genes
        """
        df.loc[:, "chromosome"] = df["Genomic Location"].str.split(":").str[0]
        df.loc[:, "chromosome"] = df["chromosome"].str.split("chr").str[1]
        df.loc[:, "Genomic Location"] = df["Genomic Location"].str.split(":").str[1]
        df.loc[:, "start"] = df["Genomic Location"].str.split("-").str[0]
        df.loc[:, "end"] = df["Genomic Location"].str.split("-").str[1]

        df = df.drop("Genomic Location", axis=1)

        return df

    def combine_PMID(self, df: pd.DataFrame) -> pd.DataFrame:

        df["PMID list Triplosensitivity"] = df[
            [
                "Triplosensitivity PMID1",
                "Triplosensitivity PMID2",
                "Triplosensitivity PMID3",
                "Triplosensitivity PMID4",
                "Triplosensitivity PMID5",
                "Triplosensitivity PMID6",
            ]
        ].values.tolist()

        df["PMID list Haploinsufficiency"] = df[
            [
                "Haploinsufficiency PMID1",
                "Haploinsufficiency PMID2",
                "Haploinsufficiency PMID3",
                "Haploinsufficiency PMID4",
                "Haploinsufficiency PMID5",
                "Haploinsufficiency PMID6",
            ]
        ].values.tolist()

        # Remove Nan values from lists
        # https://stackoverflow.com/questions/58664742/how-to-removed-nan-values-from-a-list-build-with-row-values-in-pandas-dataframe
        df["PMID list Haploinsufficiency"] = df["PMID list Haploinsufficiency"].apply(
            lambda column: [value for value in column if not pd.isna(value)]
        )
        df["PMID list Triplosensitivity"] = df["PMID list Triplosensitivity"].apply(
            lambda column: [value for value in column if not pd.isna(value)]
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """
        df["color"] = COLOR["Triplosensitivity&Haploinsufficiency"]

        df.loc[
            (
                (
                    (
                        df["Haploinsufficiency Description"]
                        != "Dosage sensitivity unlikely"
                    )
                    & (df["Haploinsufficiency Description"] != "Not yet evaluated")
                    & (df["Haploinsufficiency Description"] != "No evidence available")
                )
                & (
                    (
                        df["Triplosensitivity Description"]
                        == "Dosage sensitivity unlikely"
                    )
                    | (df["Triplosensitivity Description"] == "Not yet evaluated")
                    | (df["Triplosensitivity Description"] == "No evidence available")
                )
            ),
            "color",
        ] = COLOR["Haploinsufficiency"]

        df.loc[
            (
                (
                    (
                        df["Triplosensitivity Description"]
                        != "Dosage sensitivity unlikely"
                    )
                    & (df["Triplosensitivity Description"] != "Not yet evaluated")
                    & (df["Triplosensitivity Description"] != "No evidence available")
                )
                & (
                    (
                        df["Haploinsufficiency Description"]
                        == "Dosage sensitivity unlikely"
                    )
                    | (df["Haploinsufficiency Description"] == "Not yet evaluated")
                    | (df["Haploinsufficiency Description"] == "No evidence available")
                )
            ),
            "color",
        ] = COLOR["Triplosensitivity"]

        df.loc[
            (
                (df["Haploinsufficiency Description"] == "Dosage sensitivity unlikely")
                & (df["Triplosensitivity Description"] == "Dosage sensitivity unlikely")
            ),
            "color",
        ] = COLOR["Dosage sensitivity unlikely"]

        df.loc[
            (
                (
                    (df["Haploinsufficiency Description"] == "Not yet evaluated")
                    | (df["Haploinsufficiency Description"] == "No evidence available")
                )
                & (
                    (df["Triplosensitivity Description"] == "Not yet evaluated")
                    | (df["Triplosensitivity Description"] == "No evidence available")
                )
            ),
            "color",
        ] = COLOR["Unknown"]

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from ClinGen, see more at: https://ftp.clinicalgenome.org/README"
            + ";"
            + "Gene Symbol: "
            + df["#Gene Symbol"]
            + ";"
            + "Haploinsufficiency Description: "
            + df["Haploinsufficiency Description"]
            + ";"
            + "PMID list Haploinsufficiency: "
            + df["PMID list Haploinsufficiency"]
            + ";"
            + "Haploinsufficiency Disease ID: "
            + df["Haploinsufficiency Disease ID"]
            + ";"
            + "Triplosensitivity Description: "
            + df["Triplosensitivity Description"]
            + ";"
            + "PMID list Triplosensitivity: "
            + df["PMID list Triplosensitivity"]
            + ";"
            + "Triplosensitivity Disease ID: "
            + df["Triplosensitivity Disease ID"]
            + ";"
            + "Date Last Evaluated: "
            + df["Date Last Evaluated"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> str:

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
