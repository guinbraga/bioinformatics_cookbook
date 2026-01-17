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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Acessing NCBI database for data retrieval

```python
from Bio import Entrez, SeqIO
Entrez.email = "guilhormo.47@gmail.com"
```

```python
handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]')
rec_list = Entrez.read(handle)
if int(rec_list['RetMax']) < int(rec_list['Count']):
    handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]',
                            retmax=rec_list['Count'])
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
