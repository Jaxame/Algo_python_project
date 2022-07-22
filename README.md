# SequencesToSPSS

## Description

Within the framework of our Algorithmics unit, we developed this tool to generate a spectrum-preserving string set (SPSS) of a k-mers set extracted from sequencing data. 
This procedure is mostly conducted greedily while providing performance data. Furthermore, the query time diffences with and without the implementation of an FM-index is highlighted.

## Requirements

The core script we developped is SequencesToSPSS.py (credits: J. Ame & C. Brottier).

It uses supplementary modules to output functions runtimes and implement the FM-index: 
- timer.py (credits: C. Lemaitre & P. Peterlongo)
- Burrows_wheeler_minimal.py (credits: C. Lemaitre & P. Peterlongo)
- tools_karkkainen_sanders.py (credits: Romain Brixtel & Julien Clément)

----

## Usage

### Parameters


| Parameters | Description | Required |
|----------|:-------------:|------:|
| -i | set of reads [fasta] | Yes |
| -g | reference genome [fasta] | Yes |
| -t | solidity (abundance) threshold [int]| No (default = 2) |
| -k | k-mer size [int] | No (default = 31) |


### Example command line


	python3 SequencesToSPSS.py -i ecoli_sample_reads.fasta -g ecoli_genome_150k.fa -t 2 -k 31   


## Output

The script displays the analytic data directly in the standard ouput, as shown below (example with k=31 and t=2).

 * Example output:
	
		OUT TIME_SELECTING_KMERS=2.23
		Nb of SOLID k-mers=157975
		OUT TIME_SPSS_CONSTRUCTION=2.43
		OUT |SPSS(K)|=174588
		OUT #SPSS(K)=663
		Compression=98.16%
		OUT TIME_FN_WITHOUT_INDEX=60.03
		OUT TIME_FN_WITH_INDEX=5.75
		#FN=3634
		OUT TIME_FP_WITHOUT_INDEX=55.91
		OUT TIME_FP_WITH_INDEX=6.05
		#FP=8219
		OUT TIME_SPSS_MAX_COMPRESSION=0.05
		OUT |min SPSS(K)|=163319
		Max compression=91.82%
		Total run time = 134.81 seconds


### Output explanation:
 
 - OUT TIME_SELECTING_KMERS : time to count the solid canonical k-mers and store them (in seconds)
 - Nb of SOLID k-mers : total number of solid k-mers
 - OUT TIME_SPSS_CONSTRUCTION : time to create the SPSS from the canonical solid k-mers (in seconds)
 - OUT |SPSS(K)| : total number of characters of the SPSS
 - OUT #SPSS(K) : number of distinct sequences concatenated into the SPSS
 - Compression : compression rate of the SPSS (% of the straight distinct sequences concatenation)
 - OUT TIME_FN_WITHOUT_INDEX : False Negatives calculation time without an index (in seconds)
 - OUT TIME_FN_WITH_INDEX : False Negatives calculation time using the FM-index (in seconds) 
 - #FN : total number of False Negatives
 - OUT TIME_FP_WITHOUT_INDEX : False Positives calculation time without index (in seconds)
 - OUT TIME_FP_WITH_INDEX : False Positives calculation time using the FM-index  (in seconds)
 - #FP : total number of False Positives
 - OUT TIME_SPSS_MAX_COMPRESSION : Maximum compression time of the SPSS
 - OUT |min SPSS(K)| : total number of characters of the minimal SPSS (not used to query FN and FP)
 - Max compression : Maximum compression rate of the SPSS
 - Total runtime : in seconds

