---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: bio_genbank
    language: python
    name: bio_genbank
---

# Accessing NCBI database for data retrieval

NCBI stands for National Center for Biotechnology Information.
They possess many databases for biological data, as well as journals
and articles of a wide variety of topics. It is common for sequencing
data to be uploaded to the website, as well as other data from 
biological experiments.

Biopython has a built-in object to access NCBI and retrieve data from 
their databases. That is the `Entrez` object, which we'll explore below:

```python
from Bio import Entrez, SeqIO
Entrez.email = "guilhormo.47@gmail.com" # setting an email is good practice
```
## `EInfo`: obtaining information about Entrez databases:

EInfo provides interesting information about each of Entrez' databases,
such as field index term counts, last update, available links, and 
databases' names. We'll obtain here all databases' names:

```python
stream = Entrez.einfo()
result = stream.read()
stream.close()
```

The `result` variable now contains a list of databases in XML format (
as is standard for Entrez search return types):

```python
result
```

We can parse this result with the `Entrez` object, returning a python object:


```python
stream = Entrez.einfo()
record = Entrez.read(stream)
```

This returns a dictionary with one key, `'DbList'`. This key value is the
list of all available databases:

```python
record['DbList']
```

For each database, we can use `einfo()` again to gain more information:

```python
stream = Entrez.einfo(db='nucleotide')
record = Entrez.read(stream)
record['DbInfo'].keys()
```

```python
record['DbInfo']['Count']
```

One really useful data from this object return from `einfo(db)` is the 'FieldList',
which gives us a list of possible search fields to use with ESearch:

```python
for field in record['DbInfo']['FieldList']:
  print(f"{field['Name']}, {field['FullName']}, {field['Description']}")
```

This is used to restrict our search within a specific field, such as with
`Jones[AUTH]` to limit our search to the author field, or `Sanger[AFFL]` to
limit out search to author's affiliated to the Sanger Centre.

```python
handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]')
rec_list = Entrez.read(handle)
if int(rec_list['RetMax']) < int(rec_list['Count']):
    handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]', retmax=rec_list['Count'])
    rec_list = Entrez.read(handle)
```

```python
id_list = rec_list['IdList']
hdl = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb', retmax=rec_list['Count'])
```

```python
recs = list(SeqIO.parse(hdl, 'gb'))
```

```python
for rec in recs:
    if rec.name == 'KM288867':
        break
print(rec.name)
print(rec.description)
```

```python
for feature in rec.features:
    if feature.type == 'gene':
        print(feature.qualifiers['gene'])
    elif feature.type == 'exon':
        loc = feature.location
        print('Exon', loc.start, loc.end, loc.strand)
    else:
        print('not processed:\n%s' % feature)
```

```python
for name, value in rec.annotations.items():
  print(f'{name}={value}')
```

```python
print(rec.seq)
```

## Performing basic sequence analysis

We will be using the human **lactase** gene as an example. We will use the
same methods as above for retrieving it.

```python
from Bio import Entrez, SeqIO, SeqRecord
Entrez.email="guilhormo.47@gmail.com"
hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='gb') # lactase gene
gb_rec = SeqIO.read(hdl, 'gb')
```

Now that we have the record, we'll extract the gene sequence, even though it
contains more than that.

```python
for feature in gb_rec.features:
  if feature.type == 'CDS':
    location = feature.location # Note translation existing
cds = SeqRecord.SeqRecord(gb_rec.seq[location.start:location.end], 'NM_002299', description='LCT CDS only')
```

Now our sequence is available in the Biopython sequence record.

We'll now proceed by saving the sequence to a FASTA file:

```python
from Bio import SeqIO
with open('example.fasta', 'w') as w_hdl:
  SeqIO.write([cds], w_hdl, 'fasta')
```

The Seq.IO write method takes a list of sequences to write. In our case it is just one,
but we should be careful - if we plan to write lots of sequences, an iterator
will allocate less memory from our computer.

In most situations we will have the sequence loaded on the disk and will
want to read it. For that, we have the `SeqIO.read()` method:

```python
recs = SeqIO.parse('example.fasta', 'fasta')
for rec in recs:
  seq = rec.seq
  print(rec.description)
  print(seq[:10])
```

As we now have pure DNA, we can transcribe it:

```python
rna = seq.transcribe()
print(rna)
```

And then, we can translate the gene into a protein:

```python
prot = seq.translate()
print(prot)
```
