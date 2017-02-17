# kmer extracting using quality

* Need first using "import module biopython" in MSU HPCC
* freq_dict is the file that I saved that contained a python dict object to chaek kmer frequency of all dataset, right now I only have 11mer one.
	
	usage: python FastaProcess.py [-h] [-freq_dict DICT] [-freq FREQ]
			       k template_fastq insertion_fastq output_fasta

	positional arguments:
	  k                length of kmer
	  template_fastq   the fastq file that contain sequence we want to extract
			   kmer
	  insertion_fastq  the fastq file that contain insertion quality
	  output_fasta     the fasta file to store the extracted kmer

	optional arguments:
	  -h, --help       show this help message and exit
	  -freq_dict DICT  path that store frequency dict
	  -freq FREQ       the threshold to filter frequency
	
*use of bowtie2
	* first need to build index on the referece file, here I used pair1_2.fasta as example (bowtie2 only take fasta, so every pair I should have at least one fasta file)
		
		bowtie2-build XXX.fasta XXX_index_name
		
	* bowtie2 command I used to generate 1 error for the kmer
		
		bowtie2 -f -N 1 -L 6 -a -p 2 --score-min L,-0.8,-0.8 --norc -x index_path  kmer_fasta_path -S sam_file_to_save_alignment
		
