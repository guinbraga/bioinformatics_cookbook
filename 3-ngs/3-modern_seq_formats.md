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

<!-- ```python -->
<!-- !rm -f SRR003265.filt.fastq.gz 2>/dev/null -->
<!-- !wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz -->
<!-- ``` -->

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
We can see that there are no errors until the 25th position.
This is due to the 1000 genomes filtering rule that enforces
no N calls until the 25th position. We can notice from the
graph, however, that the distribution of N bases seems to not
be uniform, and thus is position dependent.
What about the quality of our reads?

### Studying quality of reads (phred scores) based on position

```python
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
len_recs = 0
avg_quality = defaultdict(float)
for rec in recs:
  len_recs += 1
  phred = rec.letter_annotations['phred_quality']
  for pos in range(len(rec.seq)):
    avg_quality[pos+1] += phred[pos]
for pos in range(len_recs):
  avg_quality[pos+1] /= len_recs
```

```python
positions = range(1, 52)
fig, ax = plt.subplots(figsize=(16,9))
ax.plot(positions, [avg_quality[x] for x in positions])
fig.suptitle('Average read quality per position')
ax.set_xlim(1, 52)
ax.set_ylim(0,50)
ax.set_xlabel('Read position')
ax.set_ylabel('Average Quality')
```
#### The book's solution

The bioinformatics cookbook creates this graph in a fancier way. I'll copy it below:

```python
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
cnt_qual = defaultdict(int)
for rec in recs:
  for i, qual in enumerate(rec.letter_annotations['phred_quality']):
    if i < 25:
      continue
    cnt_qual[qual] += 1
tot = sum(cnt_qual.values())
for qual, cnt in cnt_qual.items():
  print(f"{qual}: {100.*cnt/tot:.2f} {cnt}")
```

```python
recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
qual_pos = defaultdict(list)
for rec in recs:
  for i, qual in enumerate(rec.letter_annotations['phred_quality']):
    if i < 25 or qual == 40:
      continue
    pos = i+1
    qual_pos[pos].append(qual)
vps = []
poses = list(qual_pos.keys())
poses.sort()
for pos in poses:
  vps.append(qual_pos[pos])
fig, ax = plt.subplots(figsize=(16,9))
sns.boxplot(data = vps, ax=ax)
ax.set_xticklabels([str(x) for x in range(26, max(qual_pos.keys()) + 1)])
ax.set_xlabel('Read distance')
ax.set_ylabel('PHRED score')
fig.suptitle('Distribution of PHRED scores as a function of read distance')

```

## Summarizing:

We learned:
- to iterate fasta files to avoid memory allocation of the whole file;
- the contents of the fastq file, which include relevant metadata such as read quality;
- how read quality works (that is, in a log-fashion where 10 = 90% accuracy)
- some attributes of the `SeqIO.parse()` method, such as the phred_quality and the sequence itself;
- how to visualise and extract relevant data such as phred_quality as a function of base position;
- the existance of 'N' bases, which mean unidentified nucleotides; 
- some basics of plotting, and the fact that we need to develop this area!
