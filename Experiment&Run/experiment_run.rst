Experiment & Run
=================
A description of sample-specific sequencing library, instrument and sequencing methods. An experiment references 1 project and 1 sample. Runs describe the files that belong to the previously created experiments.

Genernal Information
---------------------

**\*project accession** : ``project_accession``

* Definition: A valid project accession has 'CNP' prefix.

**\*sample accession** : ``sample_accession``

* Definition: A valid sample accession has 'CNS' prefix.

**\*experiment title** : ``experiment_title``

* Definition: Short description that will identify the dataset on public pages. A clear and concise formula for the title would be like: {methodology} of {organism}: {sample info}, e.g. RNA-Seq of mus musculus: adult female spleen.

Library Information & Sequencing
--------------------------------

**\*library name** : ``library_name``

* Definition: Short unique identifier for the sequencing library. Each library name MUST be unique!

**\*library strategy** : ``library_strategy``

* Definition: Sequencing technique intended for the library.
* Value:

  | WGA: Random sequencing of the whole genome following non-pcr amplification.
  | WGS: Random sequencing of the whole genome.
  | WXS: Random sequencing of exonic regions selected from the genome.
  | RNA-Seq: Random sequencing of whole transcriptome.
  | miRNA-Seq: Random sequencing of small miRNAs.
  | WCS: Random sequencing of a whole chromosome or other replicon isolated from a genome.
  | CLONE: Genomic clone based (hierarchical) sequencing.
  | POOLCLONE: Shotgun of pooled clones (usually BACs and Fosmids).
  | AMPLICON: Sequencing of overlapping or distinct PCR or RT-PCR products.
  | CLONEEND: Clone end (5', 3', or both) sequencing.
  | FINISHING: Sequencing intended to finish (close) gaps in existing coverage.
  | ChIP-Seq: Direct sequencing of chromatin immunoprecipitates.
  | MNase-Seq: Direct sequencing following MNase digestion.
  | DNase-Hypersensitivity: Sequencing of hypersensitive sites, or segments of open chromatin that are more readily cleaved by DNaseI.
  | Bisulfite-Seq: Sequencing following treatment of DNA with bisulfite to convert cytosine residues to uracil depending on methylation status.
  | Tn-Seq: Sequencing from transposon insertion sites.
  | EST: Single pass sequencing of cDNA templates.
  | FL-cDNA: Full-length sequencing of cDNA templates.
  | CTS: Concatenated Tag Sequencing.
  | MRE-Seq: Methylation-Sensitive Restriction Enzyme Sequencing strategy.
  | MeDIP-Seq: Methylated DNA Immunoprecipitation Sequencing strategy.
  | MBD-Seq: Direct sequencing of methylated fractions sequencing strategy.
  | Synthetic-Long-Read
  | ATAC-seq: Assay for Transposase-Accessible Chromatin (ATAC) strategy is used to study genome-wide chromatin accessibility. alternative method to DNase-seq that uses an engineered Tn5 transposase to cleave DNA and to integrate primer DNA sequences into the cleaved genomic DNA.
  | ChIA-PET: Direct sequencing of proximity-ligated chromatin immunoprecipitates.
  | FAIRE-seq: Formaldehyde Assisted Isolation of Regulatory Elements. reveals regions of open chromatin.
  | Hi-C: Chromosome Conformation Capture technique where a biotin-labeled nucleotide is incorporated at the ligation junction, enabling selective purification of chimeric DNA ligation junctions followed by deep sequencing.
  | ncRNA-Seq: Capture of other non-coding RNA types, including post-translation modification types such as snRNA (small nuclear RNA) or snoRNA (small nucleolar RNA), or expression regulation types such as siRNA (small interfering RNA) or piRNA/piwi/RNA (piwi-interacting RNA).
  | RAD-Seq
  | RIP-Seq: Direct sequencing of RNA immunoprecipitates (includes CLIP-Seq, HITS-CLIP and PAR-CLIP).
  | SELEX: Systematic Evolution of Ligands by EXponential enrichment.
  | ssRNA-seq: strand-specific RNA sequencing.
  | Targeted-Capture
  | Tethered Chromatin Conformation Capture
  | OTHER: Library strategy not listed (please include additional info in the “design description”).

**\*library source** : ``library_source``

* Definition: The library source specifies the type of source material that is being sequenced.

* Value:

  | GENOMIC: Genomic DNA (includes PCR products from genomic DNA).
  | GENOMIC SINGLE CELL
  | TRANSCRIPTOMIC: Transcription products or non genomic DNA (EST, cDNA, RT-PCR, screened libraries).
  | TRANSCRIPTOMIC SINGLE CELL
  | METAGENOMIC: Mixed material from metagenome.
  | METATRANSCRIPTOMIC: Transcription products from community targets.
  | SYNTHETIC: Synthetic DNA.
  | VIRAL RNA: Viral RNA.
  | OTHER: Other, unspecified, or unknown library source material (please include additional info in the “design description”).

**\*library selection**: ``library_selection``

* Definition: Method used to enrich the target in the sequence library preparation.

* Value:

  | RANDOM: Random selection by shearing or other method.
  | PCR: Source material was selected by designed primers.
  | RANDOM PCR: Source material was selected by randomly generated primers.
  | RT-PCR: Source material was selected by reverse transcription PCR.
  | HMPR: Hypo-methylated partial restriction digest.
  | MF: Methyl Filtrated.
  | MDA: Multiple displacement amplification.
  | MSLL: Methylation Spanning Linking Library.
  | cDNA: complementary DNA.
  | ChIP: Chromatin immunoprecipitation.
  | MNase: Micrococcal Nuclease (MNase) digestion.
  | DNase: Deoxyribonuclease (MNase) digestion.
  | Hybrid Selection: Selection by hybridization in array or solution.
  | Reduced Representation: Reproducible genomic subsets, often generated by restriction fragment size selection, containing a manageable number of loci to facilitate re-sampling.
  | Restriction Digest: DNA fractionation using restriction enzymes.
  | 5-methylcytidine antibody: Selection of methylated DNA fragments using an antibody raised against 5-methylcytosine or 5-methylcytidine (m5C).
  | MBD2 protein methyl-CpG binding domain: Enrichment by methyl-CpG binding domain.
  | CAGE: Cap-analysis gene expression.
  | RACE: Rapid Amplification of cDNA Ends.
  | size fractionation: Physical selection of size appropriate targets.
  | Padlock probes capture method: Circularized oligonucleotide probes.
  | Oligo-dT: enrichment of messenger RNA (mRNA) by hybridization to Oligo-dT.
  | repeat fractionation: Selection for less repetitive (and more gene rich) sequence through Cot filtration (CF) or other fractionation techniques based on DNA kinetics.
  | cDNA_oligo_dT
  | cDNA_randomPriming
  | Inverse rRNA: depletion of ribosomal RNA by oligo hybridization.
  | PolyA: PolyA selection or enrichment for messenger RNA (mRNA); should replace cDNA enumeration.
  | other: Other library enrichment, screening, or selection process (please include additional info in the “design description”).
  | unspecified: Library enrichment, screening, or selection is not specified (please include additional info in the “design description”).

**\*library layout** : ``library_layout``

* Value:

  | fragment/single
  | paired

**\*platform** : ``platform``
  **\*instrument model** : ``instrument_model``

 - Value:

    LS454
      | 454 GS
      | 454 GS 20
      | 454 GS FLX
      | 454 GS FLX+
      | 454 GS FLX Titanium
      | 454 GS Junior

    ABI_SOLID
      | AB 5500 Genetic Analyzer
      | AB 5500xl Genetic Analyzer
      | AB 5500xl-W Genetic Analysis System
      | AB SOLiD 3 Plus System
      | AB SOLiD 4 System
      | AB SOLiD 4hq System
      | AB SOLiD PI System
      | AB SOLiD System
      | AB SOLiD System 2.0
      | AB SOLiD System 3.0

    BGISEQ
      | BGISEQ-500
      | BGISEQ-50
      | BGISEQ-1000
      | BGISEQ-100

    Bionano
      | Saphyr

    DIPSEQ
      | DIPSEQ-T1
      | DIPSEQ-T5
      | DIPSEQ-T10

    DNBSEQ
      | DNBSEQ-G50(MGISEQ-200)
      | DNBSEQ-G400(MGISEQ-2000)
      | DNBSEQ-G400 FAST
      | DNBSEQ-T1
      | DNBSEQ-T5
      | DNBSEQ-T7
      | DNBSEQ-T10
      | DNBSEQ-T10×4
      | DNBSEQ-T20
      | DNBSEQ-T20×2

    CAPILLARY
      | AB 310 Genetic Analyzer
      | AB 3130 Genetic Analyzer
      | AB 3130xL Genetic Analyzer
      | AB 3500 Genetic Analyzer
      | AB 3500xL Genetic Analyzer
      | AB 3730 Genetic Analyzer
      | AB 3730xL Genetic Analyzer

    COMPLETE_GENOMICS
      | Complete Genomics

    HELICOS
      | Helicos HeliScope

    ILLUMINA
      | HiSeq X Five
      | HiSeq X Ten
      | Illumina Genome Analyzer
      | Illumina Genome Analyzer II
      | Illumina Genome Analyzer IIx
      | Illumina HiScanSQ
      | Illumina HiSeq 1000
      | Illumina HiSeq 1500
      | Illumina HiSeq 2000
      | Illumina HiSeq 2500
      | Illumina HiSeq 3000
      | Illumina HiSeq 4000
      | Illumina iSeq 100
      | Illumina NovaSeq 6000
      | Illumina MiniSeq
      | Illumina MiSeq
      | NextSeq 500
      | NextSeq 550

    ION_TORRENT
      | Ion Torrent PGM
      | Ion Torrent Proton
      | Ion Torrent S5 XL
      | Ion Torrent S5

    OXFORD_NANOPORE
      | GridION
      | MinION
      | PromethION

    PACBIO_SMRT
      | PacBio RS
      | PacBio RS II
      | Sequel
      | Sequel II

**design description**: ``design_description``

* Definition: Free-form description of the methods used to create the sequencing library; a brief 'materials and methods' section.

**library construction protocol**: ``library_construction_protocol``

* Definition: Describes the protocol by which the sequencing library was constructed.

**\*spot layout**: ``spot_layout``

* Definition: If technical reads (e.g. barcodes, adaptors or linkers) are included in the submitted raw sequences, a spot descriptor must be submitted to describe the position of the technical reads so that they can be removed.

**\*nominal size**: ``nominal_size``

* Definition: The average insert size for paired reads.

Run Information
---------------

**\*file type**: ``file_type``

* Value:
   | bam
   | cram
   | sff
   | fastq
   | PacBio_HDF5
   | bnx
   | Oxford_Nanopore

**\*file name**: ``file_name``

**\*file md5**: ``file_md5``

Reference Information
---------------------

When submitting **BAM** files of aligned reads, you must also specify an assembly - the reference genome that your reads were aligned against. You can identify your reference assembly by its accession from the NCBI, UCSC and Ensembl. If the assembly is not available from a public repository, you will need to submit your own (local) assembly in FASTA format along with your BAM file.

When submitting **CRAM** files, the references should be provided in the same manner as BAM references.

**reference accession**: ``reference_accession``

* Definition: This is only if you are submitting a bam/cram file aligned against an assembly - the reference genome in the public repository. Please provide the accession number (e.g. GRCh37) in the public repository.

**reference fasta**: ``reference_fasta``

* Definition: Please provide the name of the custom assembly fasta file used during alignment (e.g. Mouse.fasta).

**reference md5**: ``reference_md5``

* Definition: MD5 checksum for the custom assembly fasta file.
