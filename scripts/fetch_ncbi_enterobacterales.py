# Fetch every Enterobacterales BioSample with an antibiogram; one CSV row per test.
# Usage: python3 -I scripts/fetch_ncbi_enterobacterales.py data/raw
import csv, gzip, sys, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET

out = sys.argv[1]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
term = 'antibiogram[filter] AND "Enterobacterales"[orgn]'

def get(url, tries=5):
    for k in range(tries):
        try:
            time.sleep(0.4)  # 3 requests per second without an API key
            return urllib.request.urlopen(url, timeout=120).read()
        except Exception as e:
            print("retry", k, e, file=sys.stderr); time.sleep(5 * (k + 1))
    raise RuntimeError(url)

ids, start = [], 0
while True:
    x = ET.fromstring(get(E + "esearch.fcgi?db=biosample&retmax=5000&retstart=%d&term=%s" % (start, urllib.parse.quote(term))))
    batch = [i.text for i in x.iter("Id")]
    ids += batch
    total = int(x.findtext("Count"))
    start += 5000
    if start >= total: break
print(len(ids), "ids of", total, file=sys.stderr)

cols = ["accession", "organism", "host", "collection_date", "submission_date", "publication_date",
        "antibiotic", "phenotype", "sign", "measurement", "units", "method", "standard", "header"]
with gzip.open(out + "/ncbi_enterobacterales_ast.csv.gz", "wt", newline="") as f:
    w = csv.writer(f); w.writerow(cols)
    for i in range(0, len(ids), 200):
        root = ET.fromstring(get(E + "efetch.fcgi?db=biosample&id=" + ",".join(ids[i:i + 200])))
        for s in root.iter("BioSample"):
            attr = {a.get("harmonized_name") or a.get("attribute_name"): a.text for a in s.iter("Attribute")}
            for t in s.iter("Table"):
                if not (t.get("class") or "").startswith("Antibiogram"): continue
                header = [c.text or "" for c in t.findall("./Header/Cell")]
                for r in t.findall("./Body/Row"):
                    cells = dict(zip(header, [c.text or "" for c in r.findall("./Cell")]))
                    w.writerow([s.get("accession"), s.findtext(".//Organism/OrganismName"), attr.get("host"),
                                attr.get("collection_date"), s.get("submission_date"), s.get("publication_date"),
                                cells.get("Antibiotic"), cells.get("Resistance phenotype"), cells.get("Measurement sign"),
                                cells.get("Measurement"), cells.get("Measurement units"),
                                cells.get("Laboratory typing method"), cells.get("Testing standard"), "|".join(header)])
        print(i + 200, file=sys.stderr)
