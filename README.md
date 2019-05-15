# Requirements
	python 2.7+         https://www.python.org/download/releases/2.7/
	Biopython           https://biopython.org/wiki/Download
	psiblast            ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/
	makeblastdb         ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/
	mcl                 http://micans.org/mcl/

	The bin folder has already included all the required programs
	Please always keep the required programs in the bin folder of CasLocusAnno
	Make sure you have the permission to execute the required programs before running CasLocusAnno.py

# Usage
	Command: 
	python CasLocusAnno.pyc -i <inputFile> -o <OutFile>

	Options:
	-h, --help          Check help of CasLocusAnno.py
	-v, --version       Check the current release version
	-i, --input         Input file (Whole proteome sequences)
	-o, --out           Output file of result 
	-i, -o are required parameters.

	The code can be copiable, distributable, modifiable, and usable without any restrictions.
	Report bugs to <chuand@cefg.cn> or <fbguo@uestc.edu.cn>

# Directory structure
	CasLocusAnno        Main file
	├── bin
	├── CasLocusAnno.py Main program
	├── configure
	│   └── profile.ini Profile information
	├── profiles        Cas profile
	├── readme
	├── temp            A temporary file that stores data during the annotation procedure
	└── test            File that stores test data and their annotation result
		├── NC_015709.faa
		├── res
		├── res.anno1
		└── res.anno2

# noting
	The standalone version only accepts whole proteome sequence of a chromosome.
	All the sequences in the inputting file should in the order that they appear in chromosome.

# Output
	<OutFile>.anno1     Cas protein list in the first-round search
	<OutFile>.anno2     Cas protein list in the second-round search
	<OutFile>           Cas proteins in Cas loci and (sub)types

	An example of <OutFile>. The following annotation comes from NC_015709, which is downloaded from
	ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Zymomonas_mobilis_pomaceae_ATCC_29192_uid68445/NC_015709.faa
		Cas locus type: No Cas locus
		=====
		1001	gi|338707994|ref|YP_004662195.1|	7.32e-43	pfam09707.sr	cas2	cas2	CAS-I-E
		1002	gi|338707995|ref|YP_004662196.1|	3.02e-132	cd09719.sr	cas1	cas1	CAS-I-E
		1003	gi|338707996|ref|YP_004662197.1|	9.04e-37	pfam08798.sr	cas6e	cas6e	CAS-I-E,CAS-IV
		1004	gi|338707997|ref|YP_004662198.1|	3.58e-35	cd09645.sr	cas5	cas5	CAS-I-E
		1005	gi|338707998|ref|YP_004662199.1|	2.91e-115	pfam09344.sr	cas7	cas7	CAS-I-E
		1006	gi|338707999|ref|YP_004662200.1|	9.69e-20	cd09731.sr	cse2gr11	cse2gr11	CAS-I-E
		1007	gi|338708000|ref|YP_004662201.1|	3.8e-83	cd09729.sr	cas8e	cas8e	CAS-I-E
		1008	gi|338708001|ref|YP_004662202.1|	3.31e-64	COG1203.sr	cas3	cas3	CAS-I
		1008	gi|338708001|ref|YP_004662202.1|	2.05e-21	cd09641.sr	cas3HD	cas3HD	CAS-I
		Cas locus type: CAS-I-E
		=====
		Cas locus type: No Cas locus
		=====
		Cas locus type: No Cas locus
		=====
	We used MCCS further to filter the initial annotation,
	Cas proteins that may not constitute a Cas locus or may be false discovery were discarded.
	This is the reason why you see "Cas locus type: No Cas locus".
	There is one Cas locus belonging to CAS-I-E in the chromosome 