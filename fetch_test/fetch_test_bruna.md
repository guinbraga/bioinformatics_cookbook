---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: md,ipynb
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

# My own fetch test on bruna's data:

First, I'll initialize the Entrez object, set my email and check the
databases' FieldList for search:


```python
from Bio import Entrez
Entrez.email = 'guilhormo.47@gmail.com'
```

```python
handle = Entrez.einfo(db='bioproject')
result = Entrez.read(handle)
result.keys()
result['DbInfo'].keys()
result['DbInfo']['FieldList']
for key in result['DbInfo']['FieldList']:
  print(f"{key['Name']}: {key['FullName']}, {key['Description']}")
```

It looks like we'll search for the `PRJNA788342[PRJA]`

```python
handle = Entrez.esearch(db='bioproject', term='PRJNA788342[PRJA]')
result = Entrez.read(handle)
project_id = result['IdList'][0]
result
```

```python
# Remember to close the handle!
handle.close()
```

Now that we've found the project ID in the bioproject database, we'll
search for it's data in the SRA database

```python
handle = Entrez.elink(dbfrom='bioproject', db='sra', id=project_id)
result = Entrez.read(handle)
result
result[0].keys()
```

```python
result[0]['LinkSetDb'][0].keys()
```

```python
result[0]['LinkSetDb'][0]['Link']
```


```python
handle.close()
sra_ids = [id['Id'] for id in result[0]['LinkSetDb'][0]['Link']]
len(sra_ids)
with open('sra_ids.txt', 'w') as output:
  for id in sra_ids:
    output.write(f'{id}\n')
```

The above object contains the ID for each sample taken during the experiment.
Each sample contained thousands of reads, but the authors rarefied them to 
about 50 thousand reads per sample.
Next up (probably tomorrow) we'll try and fetch some of those FASTQ files!
