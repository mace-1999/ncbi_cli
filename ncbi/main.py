from ncbi.ictv_scrape import GetVirusInfo
import click
from ncbi.abstract_get import abs


@click.group()
@click.version_option(message="%(prog)s version %(version)s")
def cli():
    pass


@cli.command()
@click.option("--search", required=True, help="Search term for returned abstracts")
@click.option("--email", required=True, help="Your email address for the NCBI API")
@click.option("--bs", default=5, help="Download abstracts in batches of n", type=int, show_default=True)
@click.option("--total_size", default=20, help="Total n abstracts to download", type=int, show_default=True)
def abstract(search, email, bs, total_size):
    abs(search, email, bs, total_size)


@cli.command()
@click.option("--name", required=True, help="Virus name to generate report")
@click.option("--email", required=True, help="Your email address for the NCBI API", type=str)
@click.option("--csv", default="F", help="Enter T / F to keep raw csv file", show_default=True, type=str)
def virus(name, email, csv):
    GetVirusInfo(name,email=email,  csv=csv)


if __name__ == "__main__":
    cli()
