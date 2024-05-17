import sys

import requests
from bs4 import BeautifulSoup
from requests.exceptions import HTTPError
import re
import pandas as pd
from Bio import SeqIO, Entrez
import os
import plotly.express as px
from fpdf import FPDF
import click
from urllib.error import URLError

class GetVirusInfo:
    def __init__(self, virus: "str", csv, email: "str"):
        self.virus = virus
        self.csv = csv
        self.filename = f"{self.virus}.gb"
        self.plot_filename = f"{self.virus}_plot.png"
        self.df_filename = f"{self.virus}.csv"
        self.pdf_out = f"{self.virus}_report.pdf"
        self.WIDTH = 210
        self.HEIGHT = 297
        self.soup = None
        self.email = email
        self.download_genbank()
        self.dict = {"Name": [], "Accession": [], "Length": []}
        self.name()
        self.get_accession()
        self.length()
        self.df = pd.DataFrame(self.dict).sort_values(by="Length", ascending=False)
        self.get_figure()
        self.checked_scrape = self.get_soup()
        if self.soup:
            self.summary = self.get_summary_art()
            self.tree = self.get_fam_tree()
            self.create_report()
            self.clean_dir()
            click.launch(self.pdf_out)
        else:
            click.echo(click.style("Unable to generate report", fg='red'))
            click.echo(click.style("Generating Raw files...", fg='red'))
            self.clean_dir(saved=True)

    def get_soup(self):
        ictv_url = f"https://ictv.global/report/chapter/{self.virus}/{self.virus}"

        try:
            response = requests.get(ictv_url)

            # If the response was successful, no Exception will be raised
            response.raise_for_status()
            response_text = response.text
        except HTTPError as http_err:
            click.echo(click.style(f'HTTP error occurred: {http_err}', fg='red'))
            return False

        else:
            # Successful response results in parse through html.parser
            print("Success")
            self.soup = BeautifulSoup(response_text, "html.parser")
            return True

    def get_summary_art(self):
        """
                    A function to parse viruses and get their summary from the ictv website

                    :param soup: the name of a virus family
                    :return: string of virus summary

        """
        summary_art = self.soup.find_all(lambda tag: tag.name == 'p' or tag.name == "h2")
        summary_art_list = [i.getText().strip() for i in summary_art]
        # print(news_arts_list)
        num = summary_art_list.index("Summary")
        summary_article_fin = summary_art_list[num + 1]
        summary_article_fin = re.sub(" [(\[].*?[)\]]", "", summary_article_fin).strip()
        summary_article_fin = summary_article_fin.encode('latin-1', 'replace').decode('latin-1')
        return summary_article_fin

    # soup = get_soup("paramyxoviridae")
    # get the family tree

    def get_fam_tree(self):
        menu = self.soup.find_all(class_="menu-item menu-item--expanded")
        menu_list = []
        for i in menu:
            menu_items = i.find_all(hreflang="en")
            menu_list += [i.getText().strip() for i in menu_items]
        return menu_list

    def download_genbank(self):
        Entrez.email = self.email  # Always tell NCBI who you are
        try:
            handle = Entrez.egquery(term=self.virus + "[ORGN] AND refseq[FILT] AND biomol_genomic[PROP]")
            record = Entrez.read(handle)
            max_hits = 50
            for row in record["eGQueryResult"]:
                if row["DbName"] == "nuccore":
                    max_hits = int(row["Count"])

                    print('Total Sequences... ')

                    total_sequences = row["Count"]
                    print(total_sequences)

        except URLError:
            click.echo(click.style("Egquery down - defaulting max hits to 200.", fg='red'))
            max_hits = 200


        search_handle = Entrez.esearch(db='nuccore',
                                       term=self.virus + "[ORGN] AND refseq[FILT] AND biomol_genomic[PROP]",
                                       usehistory='y',
                                       retmax=max_hits)

        search_results = Entrez.read(search_handle)
        search_handle.close()

        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        gi_list = search_results['IdList']
        print('ID List Length... ')
        print(len(gi_list))
        # print('Total Sequences... ')
        # print(search_results['Count'])

        print('Downloading Results... ')

        fetch_handle = Entrez.efetch(db='nuccore',
                                     rettype="gb",
                                     retmode="text",
                                     retmax=max_hits,
                                     webenv=webenv,
                                     query_key=query_key)

        data = fetch_handle.read()
        fetch_handle.close()

        out_handle = open(self.filename, 'w')
        try:
            out_handle.write(data)
        except TypeError:
            click.echo(click.style(f'Error parsing data - check if virus family exists.', fg='red'))
            sys.exit()

        out_handle.close()

        click.echo(click.style("Download Complete!", fg='green'))

    def name(self):
        micro = open(self.filename, 'r')
        for seq_record in SeqIO.parse(micro, 'genbank'):
            self.dict["Name"] += [seq_record.annotations['organism']]

    def get_accession(self):
        micro = open(self.filename, 'r')
        for seq_record in SeqIO.parse(micro, 'genbank'):
            self.dict["Accession"] += [seq_record.annotations['accessions'][0]]

    def length(self):
        micro = open(self.filename, 'r')
        for seq_record in SeqIO.parse(micro, 'genbank'):
            self.dict["Length"] += [len(seq_record)]

    def get_figure(self):

        fig = px.bar(self.df, x='Accession', y='Length', color="Length", text="Length")
        fig.update_traces(texttemplate='%{text:.2s}')
        fig.update_layout(width=550, height=500,
                          title_text=f"Length range of {self.virus.title()} virus genomes (bases)",
                          font_color="black")
        fig.update_xaxes(showticklabels=False)
        fig.write_image(self.plot_filename)

    def create_report(self):

        pdf = FPDF('P', 'mm', 'A4')
        pdf.add_page()
        pdf.set_font('Arial', 'B', 16)
        pdf.ln(h="50")
        pdf.ln(h="50")
        pdf.ln(h="50")
        pdf.cell(self.WIDTH - 30, 10, f'{self.virus.title()}', align="C")
        pdf.ln(h="50")
        pdf.cell(self.WIDTH - 30, 10, 'Summary', align="L")
        pdf.ln(h="50")
        pdf.set_font("Arial", "", 10)
        pdf.write(5, self.summary)

        pdf.ln(h="")
        pdf.ln(h="")
        pdf.set_font('Arial', 'B', 16)
        pdf.cell(self.WIDTH - 30, 10, 'Family Tree', align="L")
        pdf.ln(h="")
        pdf.set_font("Arial", "B", 10)
        for x, word in enumerate(self.tree):
            # print(x)
            pdf.write(5, word)
            pdf.write(5, " ")
            try:
                if "Subfamily" in self.tree[x + 1]:
                    pdf.ln()
            except IndexError:
                continue

        pdf.image(self.plot_filename,
                  y=120,
                  x=0, w=self.WIDTH - 10,
                  h= self.HEIGHT - 100)

        pdf.output(self.pdf_out, 'F')

    def clean_dir(self, saved=False):

        os.remove(self.filename)
        if saved:
            self.df.to_csv(self.df_filename)
            click.echo(f"Files saved as {self.df_filename} & "
                       f"{self.plot_filename}")
        elif self.csv == "T":
            self.df.to_csv(self.df_filename, index=False)
            click.echo(f"CSV saved as {self.df_filename}")
            os.remove(self.plot_filename)
        else:
            os.remove(self.plot_filename)
            print("Cleaning directory...")
