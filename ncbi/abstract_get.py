from Bio import Entrez
import click


def abs(search, email, bs, total_size):
    Entrez.email = email  # Always tell NCBI who you are
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed", term=search, usehistory="y", sort="relevance"
        )
    )
    count = int(search_results["Count"])
    click.echo(f"Found {count} results for {search}")
    if total_size >= count:
        click.echo(f"Downloading {count} Results...")
    else:
        click.echo(f"Downloading {total_size} Results...")
    # click.pause("Press any key to download ?")
    with click.progressbar(length=total_size, empty_char="-",
                           fill_char="*", bar_template='[%(bar)s]  %(info)s') as bar:
        out_handle = open(f"{search}_abstracts.txt", "w")
        for start in range(0, count, bs):
            end = min(count, start + bs)
            if start + 1 >= total_size:
                click.echo(click.style("Max size reached", fg='red'))
                bar.update(total_size)
                break

            # print(f"Going to download record {start + 1} to {end}")
            fetch_handle = Entrez.efetch(
                db="pubmed",
                rettype="abstract",
                retmode="text",
                retstart=start,
                retmax=bs,
                webenv=search_results["WebEnv"],
                query_key=search_results["QueryKey"],
            )
            data = fetch_handle.read()
            fetch_handle.close()
            out_handle.write(data)
            if total_size >= count:
                bar.update(total_size)
            else:
                bar.update(end)
    out_handle.close()

    click.echo(click.style("Download Complete!", fg='green'))
    click.launch(f"{search}_abstracts.txt")
