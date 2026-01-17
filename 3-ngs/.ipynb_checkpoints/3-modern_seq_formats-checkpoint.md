---
jupyter:
  jupytext:
    cell_metadata_filter: -all
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

# Working with modern sequence formats

In this notebook, we'll use the Humans 1,000 genomes project data to explore
the usage of FASTQ files with their quality score for each base.

One of the main challenges of using NGS data is the raw size of the data.
Due to it's large nature, we must be aware of ways to handle disk space and
backup policies.

```python
!rm -f SRR003265.filt.fastq.gz 2>/dev/null
!wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz
```

```python
import gzip
from Bio import SeqIO
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
rec = next(recs)
print(rec.id, rec.description, rec.seq)
print(rec.letter_annotations)
```

note that here `recs` holds an iterator with the contents of the FASTQ file. This means
that if we were to convert this iterator to a list, it might eat our whole RAM. **The 
safest way to manipulate data in a FASTQ file is to either do all manipulations in a 
single iteration or open and close the file multiple times over.**

### taking a look at the distribution of nucleotide reads:



```python
from collections import defaultdict
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
cnt = defaultdict(int)
for rec in recs:
  for letter in rec.seq:
    cnt[letter] += 1
tot = sum(cnt.values())
for letter, count in cnt.items():
  print(f"{letter}: {(100.*count/tot):.2f} {count}")
```

Here, N's represent faulty data, or base pairs that could not be identified.

### plotting the distribution of unidentified base pairs according to their reading position:

before doing this recipe, I wondered: whats the length of each sequence in the record?


```python
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
avg_len = 0
seq_cnt = 0
for rec in recs:
  avg_len += len(rec.seq) 
  seq_cnt += 1
avg_len /= seq_cnt
print(avg_len)
```
For some reason, these are really small sequences. I could not find in the book the reason why,
nor what gene is this sequence from, except it is from Youruban people.

Let's continue with the plot of the distance until the unidentified bases:

```python
import seaborn as sns
import matplotlib.pyplot as plt
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
n_cnt = defaultdict(int)
for rec in recs:
  for i, letter in enumerate(rec.seq):
    if letter == 'N':
      n_cnt[i+1] += 1
seq_len = max(n_cnt.keys())
positions = range(1, seq_len + 1) # these will serve as x axis
fig, ax = plt.subplots(figsize=(16, 9))
ax.plot(positions, [n_cnt[x] for x in positions])
fig.suptitle('Number of N calls as a function of the distance from the start of the sequencer read')
ax.set_xlim(1, seq_len)
ax.set_xlabel('Read distance')
ax.set_ylabel('Number of N Calls')
```
