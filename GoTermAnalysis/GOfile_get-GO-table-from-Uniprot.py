from requests.adapters import HTTPAdapter, Retry
import re
import requests
import pandas as pd

def download_organism_go_terms(tax_id, outfile):
    re_next_link = re.compile(r'<(.+)>; rel="next"')

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    #url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cgo&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
    url = f"https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cgene_oln%2Corganism_name%2Clength%2Cgo_id%2Cft_transmem%2Cft_intramem&format=tsv&query=%28%28taxonomy_id%3A{tax_id}08%29%29&size=500"
    progress = 0
    print("Starting to Download GO Terms")
    with open(outfile, "w") as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'Downloaded {progress} / {total}')
    return 0


def generateTagtoGOMapping(mastertable):
    #uniprotgo = "PATH/TO/TABLE"
    uniprotgo = mastertable
    df = pd.read_csv(uniprotgo, sep="\t")

    df = df.rename({"Gene Names": "locus_tag"}, axis=1)
    df = df[~df["locus_tag"].isna()]
    #df = df.rename({"Gene Names (ordered locus)": "locus_tag"}, axis=1)
    #df = df[~df["locus_tag"].isna()]
    symbols = df[["Entry", "locus_tag"]]
    df["GOTerm"] = df["Gene Ontology IDs"].str.split('; | |/')
    df = df[["locus_tag", "GOTerm"]]
    df = df.explode("locus_tag").explode("GOTerm")
    df = df[~df["GOTerm"].isna()]
    df.to_csv("output_go_terms.log", sep="\t", index=False) # Output durch pfad zur neuen datei ersetzen
    symbols.to_csv("output_symbols.log", sep="\t", index=False)

def generateTagtoGOMapping_withUTR(mastertable):
    """
    same as generateTagtoGOMapping but adds a 5'UTR and a 3'UTR to each gene
    """

    def addUTRs(df):
        """
        for each gene (row in the dataframe column called locus_tag), add one row for 3'-UTR and one for 5'UTR of the same gene
        the locus tags becomes "3'UTR-{locus_tag}" or "5'UTR-{locus_tag}", all the other entries of the row stay the same
        """
        UTRs = ["3'UTR-", "5'UTR-"] # what gets added, in this cases the 3' and 5' UTRs but can be replaced by anything
        for index, row in df.iterrows():
            for utr in UTRs:
                new_row = {}
                for column_name in df:
                    if column_name == "locus_tag":
                        new_row[column_name] = utr + row[column_name]
                    else:
                        new_row[column_name] = row[column_name]
                df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
        return df

    #uniprotgo = "PATH/TO/TABLE"
    uniprotgo = mastertable
    df = pd.read_csv(uniprotgo, sep="\t")
    df = df.rename({"Gene Names": "locus_tag"}, axis=1)
    df = df[~df["locus_tag"].isna()]
    #df = df.rename({"Gene Names (ordered locus)": "locus_tag"}, axis=1)
    #df = df[~df["locus_tag"].isna()]
    symbols = df[["Entry", "locus_tag"]]
    df["GOTerm"] = df["Gene Ontology IDs"].str.split('; | |/')
    df = df[["locus_tag", "GOTerm"]]
    df = addUTRs(df)
    symbols = addUTRs(symbols)
    df = df.explode("locus_tag").explode("GOTerm")
    df = df[~df["GOTerm"].isna()]
    df.to_csv("output_go_terms.log", sep="\t", index=False) # Output durch pfad zur neuen datei ersetzen
    symbols.to_csv("output_symbols.log", sep="\t", index=False)

##############################################################################################################

outfile = "synechocystisGOterms.txt"
download_organism_go_terms(11117, outfile)
generateTagtoGOMapping_withUTR(outfile)