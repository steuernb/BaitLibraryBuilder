# BaitLibraryBuilder

BaitLibraryBuilder will create baits from source sequences. Whenever a bait is created redundant sub-sequences in all source sequences are soft-masked. Soft-masked areas are ignored for further bait design except to elongate a non-masked sequence up to bait length. The criterium for elongation into soft-masked sequence is set by -f, which is 1.5 by dedfault. As long as the non-masked sequence is longer than 1.5th fraction (75%) of bait length.



## Data preparation

### Prepare input sequences
The source sequences for generating baits have to be collected. 
Sequences within that should not be part of baits, such as repeats have to be replaced with stretches of "N".

If sequences from different sources are to be treated differently while designing the baits, a prefix should be added to the identifyer followed by a "_".

Sequences from all sources have to be in fasta format and in a single file.

### Blast
NCBI blast is needed for this step.

One input file for the bait design is a blast of all source sequences against themselves. This is to remove redundancy. The output format for the blast should be xml. This can be set in NCBI blast with `-outfmt 5`.

### Identity and input order
If different sources, for example different species, are used in the design process, the removal of redundant sequences can be performed with adjusted identity threshold. For example, source sequences from the main species can be included and redundancy is only removed if it is with very high identity (99%) whereas less important sources are only included if they are very different from already included data (such as 80%). The identity thresholds can be set in a TSV file where the first column is the prefix of a set of source sequences and the second column is the identity percentage. The order of prefixes in the TSV file also determines the order in which the different input source sets are processed.

## Running BaitLibraryBuilder

Java runtime environment is needed to run the BaitLibraryBuilder

```
java -jar BaitLibraryBuilder.jar -a baitlength -b blast.xml -d threshold.tsv -i concatednatedInput.fasta -o baits.fasta -v overlap
```

Parameter | Argument type | default value | Description
--- | --- | ---- | ---
-i | FASTA | mandatory |  A file in fasta format containing the source sequences
-o | File |mandatory | The output file that will be generated. Format is fasta
-b | XML |mandatory | The output of the self blast of the sources. Format is XML (`-outfmt 5` in blast)
-a | int| 120 |  length of the baits
-d | File / double | 95 | Either a double value denoting the identity thresold for removal of redundant sequences or a TSV file with identity threshold for each source data set
-v | int | 0 | number of base pairs that baits can overlap. If this is negative there is a gap between baits
-f | float | 1.5 |  minimum fraction of bait length to ignore
-s | int | 100 | length threshold of an alignment (blast hsp) for redundancy

