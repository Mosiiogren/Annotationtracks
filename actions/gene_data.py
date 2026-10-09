import re
import io
import datetime
import requests
import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path
from st2common.runners.base_action import Action

COLOR = {
    "mane plus clinical": "144, 238, 144",
    "mane select": "120. 200, 120",
    "gencode_primary": "30, 70, 140",
    "ensembl_canonical": "30, 70, 140",
    "gencode_basic": "70, 120, 190",
    "Unknown": "153, 153, 153",
}


class GeneData(Action):
    """
    Class for selecting MANE status or Ensemble tagged transcripts/genes and their corresponding exons
    """

    def run(
        self,
        MANEfileurl: str,
        genefileurl: str,
        columns: list[str],
        attributes: list[str],
        outputfilegene: str,
        outputfileexon: str,
        outputfoldertracks: str,
        filename: str,
    ) -> tuple[bool, str]:

        succeded, results = self.get_MANE_data(MANEfileurl)
        if not succeded:
            return (False, results)
        df_MANE_data = self.filter_MANE_data(results)

        succeded, results = self.get_gene_data(genefileurl, columns, attributes)
        if not succeded:
            return (False, results)
        df_ENSEMBL_data = self.filter_gene_data(results)

        df_ENSEMBL_genes = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "gene")]
        df_transcripts = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "transcript")]
        df_exons = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "exon")]

        df_genes = self.create_MANE_genes_dataframe(
            df_transcripts, df_MANE_data, df_ENSEMBL_genes
        )

        df_exons = self.crete_MANE_exons_dataframe(df_exons, df_genes)

        df_genes.to_json(outputfilegene, orient="records")
        df_exons.to_json(outputfileexon, orient="records")

        df_genes = self.add_comments(df_genes)
        df_genes = self.add_color(df_genes)
        filename = self.create_annotationtrack_files(
            df_genes,
            outputfoldertracks,
            filename,
        )

        return (True, filename)

    def get_MANE_data(self, url: str) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving MANE data from request
        """
        succeded, results = self.get_response(url)
        if not succeded:
            return False, results
        df = pd.read_csv(
            io.StringIO(results.decode()),
            sep="\t",
            na_values="\\N",
            dtype=str,
        )

        return True, df

    def get_gene_data(
        self, url: str, columns: list[str], attributes: list[str]
    ) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving gene data from request
        """

        succeded, results = self.get_response(url)
        if not succeded:
            return False, results
        df = pd.read_csv(
            io.StringIO(results.decode()),
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
        )
        try:
            df.columns = columns
        except:
            return (False, f"Invalid columns for File {url}")

        succeded, results = self.get_attributes(
            df=df,
            attributes=attributes,
            last_column="attribute",
        )

        if not succeded:
            return False, results

        return True, results

    def get_response(self, url: str) -> tuple[bool, str]:
        """
        Function for retrieving the response from a request
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

        return (True, response)

    def get_attributes(
        self, df: pd.DataFrame, attributes: list, last_column
    ) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for extracting key value data
        """
        try:
            for attribute in attributes:
                df[attribute] = df[last_column].apply(
                    lambda x: (
                        re.findall(rf'{attribute} "([^"]*)"', x)[0]
                        if rf'{attribute} "' in x
                        else np.nan
                    )
                )
        except:
            return (False, f"Invalid attribute list {attributes}")
        df = df.drop(last_column, axis=1)

        return True, df

    def filter_MANE_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the data into the correct format
        """

        # Filter the chromosome column -> Looks like: NC_000019.10 (We want 19)
        df = df[df["GRCh38_chr"].str.contains("NC") == True].copy()
        df.loc[:, "GRCh38_chr"] = df["GRCh38_chr"].str.split(".").str[0]
        df.loc[:, "GRCh38_chr"] = df["GRCh38_chr"].str.split("_").str[1]
        df["GRCh38_chr"] = df["GRCh38_chr"].astype("int").astype("str")
        df.loc[df["GRCh38_chr"] == "23", "GRCh38_chr"] = "X"
        df.loc[df["GRCh38_chr"] == "24", "GRCh38_chr"] = "Y"

        # Remove the version from Ensemble ID and Ensembl_nuc
        df.loc[:, "Ensembl_Gene"] = df["Ensembl_Gene"].str.split(".").str[0]
        df.loc[:, "Ensembl_nuc"] = df["Ensembl_nuc"].str.split(".").str[0]
        df = df.rename(
            columns={
                "Ensembl_nuc": "transcript_id",
                "Ensembl_Gene": "gene_id",
                "chr_start": "start",
                "chr_end": "end",
                "GRCh38_chr": "chromosome",
                "symbol": "gene_name",
                "chr_strand": "strand",
            }
        )

        df = df.drop(
            df.columns.difference(
                [
                    "gene_id",
                    "gene_name",
                    "transcript_id",
                    "MANE_status",
                    "chromosome",
                    "start",
                    "end",
                    "strand",
                ]
            ),
            axis=1,
        )

        return df

    def filter_gene_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the data into the correct format
        """

        # Remove version from chromosome
        df.loc[:, "seqname"] = df["seqname"].str.split(".").str[0]
        df.rename(columns={"seqname": "chromosome"}, inplace=True)

        df = df[df["chromosome"].str.contains("Un|EBV|random|M|KI|GL") == False]
        df = df.drop(["score", "frame", "source"], axis=1)

        return df

    def create_MANE_genes_dataframe(
        self,
        df_transcripts: pd.DataFrame,
        df_MANE: pd.DataFrame,
        df_genes: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Function for slecting MANE status and Ensemble tagged transcripts
        """

        # Retrieve all transcripts of genes not included in the MANE datafile
        df_transcripts_Ensemble = df_transcripts[
            ~df_transcripts.gene_id.isin(df_MANE.gene_id)
        ]
        df_transcripts_Ensemble.loc[
            (
                (df_transcripts_Ensemble["tag"] == "Ensembl_canonical")
                | (df_transcripts_Ensemble["tag"] == "gencode_primary")
            ),
            "MANE_status",
        ] = df_transcripts_Ensemble["tag"]

        df_transcripts_Ensemble = df_transcripts_Ensemble[
            (df_transcripts_Ensemble["MANE_status"].notna())
        ]

        # Add all the transcripts together
        df_all_transcripts = pd.concat([df_transcripts_Ensemble, df_MANE])

        # Retrieve all genes not included in the transcript datafile
        df_genes_Ensemble = df_genes[~df_genes.gene_id.isin(df_all_transcripts.gene_id)]
        df_all_genes = pd.concat([df_genes_Ensemble, df_all_transcripts])

        df_all_genes.loc[df_all_genes["MANE_status"].isna(), "MANE_status"] = "-"
        df_all_genes.loc[df_all_genes["transcript_id"].isna(), "transcript_id"] = "-"

        return df_all_genes

    def crete_MANE_exons_dataframe(self, df_exons, df_genes) -> pd.DataFrame:
        """
        Function for slecting the exons belonging to the MANE status and Ensemble tagged transcripts
        """

        df_genes = df_genes[df_genes["MANE_status"].notna()]
        df_exons = df_exons[df_exons.transcript_id.isin(df_genes.transcript_id)]

        return df_exons

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Data retrieved from NCBI, see more at: https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/README and  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/README.txt"
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
            + "MANE Status: "
            + df["MANE_status"]
            + ";"
            + "Transcript ID: "
            + df["transcript_id"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for status in df["MANE_status"].unique():
            if status.lower() in COLOR.keys():
                df.loc[df["MANE_status"] == status, "color"] = COLOR[status.lower()]
            else:
                df.loc[df["MANE_status"] == status, "color"] = COLOR["Unknown"]

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
                    "color",
                ]
            ),
            axis=1,
        )
        df["start"] = df["start"].astype("float").astype("int")
        df["end"] = df["end"].astype("float").astype("int")

        df.to_csv(
            (outputfolder + filename),
            sep="\t",
            index=False,
        )

        return [filename]
