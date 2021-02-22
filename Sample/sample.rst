Sample
========

Description of biological source material; each physically unique specimen should be registered as a single sample with a unique set of attributes.

Sample classification
---------------------

non-sample terms
~~~~~~~~~~~~~~~~

**estimated size** : ``estimated_size``

* Definition: The estimated size of the genome (in bp) prior to sequencing. Of particular importance in the sequencing of (eukaryotic) genome which could remain in draft form for a long or unspecified period.
* Field Format: restricted text
* Expected value: number of base pairs
* Value syntax: {integer} bp
* Example: 300000 bp

**experimental factor** : ``experimental_factor``

* Definition: Experimental factors are essentially the variable aspects of an experiment design which can be used to describe an experiment, or set of experiments, in an increasingly detailed manner. This field accepts ontology terms from Experimental Factor Ontology (EFO) and/or Ontology for Biomedical Investigations (OBI). For a browser of EFO (v 2.95) terms, please see http://purl.bioontology.org/ontology/EFO; for a browser of OBI (v 2018-02-12) terms please see http://purl.bioontology.org/ontology/OBI
* Field Format: free text
* Expected value: text or EFO and/or OBI
* Value syntax: {termLabel} {[termID]}|{text}
* Example: time series design [EFO:EFO_0001779]

**extrachromosomal elements** : ``extrachrom_elements``

* Definition: Do plasmids exist of significant phenotypic consequence (e.g. ones that determine virulence or antibiotic resistance). Megaplasmids? Other plasmids (borrelia has 15+ plasmids)
* Field Format: restricted text
* Expected value: number of extrachromosmal elements
* Value syntax: {integer}

**investigation type** : ``investigation_type``

* Definition: Nucleic Acid Sequence Report is the root element of all MIxS compliant reports as standardized by Genomic Standards Consortium.
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'bacteria_archaea', 'eukaryote', 'metagenome', 'metagenome-assembled genome', 'metatranscriptome', 'mimarks-specimen', 'mimarks-survey', 'organelle', 'plasmid', 'single amplified genome', 'uncultivated virus genome', 'virus']

**number of replicons** : ``num_replicons``

* Definition: Reports the number of replicons in a nuclear genome of eukaryotes, in the genome of a bacterium or archaea or the number of segments in a segmented virus. Always applied to the haploid chromosome count of a eukaryote.
* Field Format: restricted text
* Expected value: for eukaryotes and bacteria: chromosomes (haploid count); for viruses: segments
* Value syntax: {integer}

**ploidy** : ``ploidy``

* Definition: The ploidy level of the genome (e.g. allopolyploid, haploid, diploid, triploid, tetraploid). It has implications for the downstream study of duplicated gene and regions of the genomes (and perhaps for difficulties in assembly). For terms, please select terms listed under class ploidy (PATO:001374) of Phenotypic Quality Ontology (PATO), and for a browser of PATO (v1.269) please refer to http://purl.bioontology.org/ontology/PATO
* Field Format: free text
* Expected value: PATO
* Value syntax: {term}
* Example: allopolyploid, polyploid

**pooling of DNA extracts** : ``pool_dna_extracts``

* Definition: were multiple DNA extractions mixed? how many?
* Field Format: free text
* Expected value: pooling status;measurement value
* Value syntax: {boolean};{float} {unit}
* Preferred unit: gram, milliliter, microliter

**project accession** : ``project_accession``

* Definition: A valid project accession has 'CNP' prefix.
* Field Format: restricted text
* Expected value: CNSA project accession

**reference for biomaterial** : ``ref_biomaterial``

* Definition: primary publication if isolated before genome publication; otherwise, primary genome report
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**sample volume or weight for DNA extraction** : ``samp_vol_we_dna_ext``

* Definition: volume (mL) or weight (g) of sample processed for DNA extraction
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millliter, gram, milligram

**study completion status** : ``study_complt_stat``

* Definition: specification of study completion status, if no the reason should be specified
* Field Format: text choice
* Expected value: study completion status
* Value syntax: {boolean};[adverse event\|non-compliance\|lost to follow up\|other-specify]
* Example: No - adverse event; No - lost to follow up; No - non-compliance; No - other; Yes.

collection event information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**altitude** : ``altitude``

* Definition: the altitude of the sample is the vertical distance between Earth's surface above sea level and the sampled position in the air.
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} m
* Preferred unit: meter

**collection date** : ``collection_date``

* Definition: The time of sampling, either as an instance (single point in time) or interval. Date/time ranges are supported by providing two dates from among the supported value formats, delimited by a forward-slash character,e.g., 2017/2019; In case no exact time is available, the date/time can be right truncated i.e. all of these are valid ISO8601 compliant times: 2008-01-23T19:23:10+00:00; 2008-01-23T19:23:10; 2008-01-23; 2008-01; 2008.
* Field Format: restricted text
* Expected value: date and time
* Value syntax: {timestamp}

**depth** : ``depth``

* Definition: depth is defined as the vertical distance below local surface, e.g. for sediment or soil samples depth is measured from sediment or soil surface, respectivly. Depth can be reported as an interval for subsurface samples
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} m
* Preferred unit: m

**elevation** : ``elev``

* Definition: the elevation of the sampling site as measured by the vertical distance from mean sea level
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: meter

**environment (biome)** : ``env_biome``

* Definition:
    | The environmental biome level are the major classes of ecologically similar communities of plants, animals, and other organisms. Biomes are defined based on factors such as plant structures, leaf types, plant spacing, and other factors like climate. Biome should be treated as the descriptor of the broad ecological context of a sample. Examples include: desert, taiga, deciduous woodland, or coral reef. EnvO (v1.53) terms listed under environmental biome can be found from the link: http://www.environmentontology.org/Browse-EnvO or https://www.ebi.ac.uk/ols/ontologies/envo
    | Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome [ENVO:00000428]. Multiple terms can be separated by one or more pipes e.g.: mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]
* Field Format: free text
* Expected value: EnvO
* Value syntax: {term}
* Example: mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]

**environment (feature)** : ``env_feature``

* Definition:
    | Environmental feature level includes geographic environmental features. Compared to biome, feature is a descriptor of the more local environment. Examples include: harbor, cliff, or lake. EnvO (v1.53) terms listed under environmental feature can be found from the link: http://www.environmentontology.org/Browse-EnvO or https://www.ebi.ac.uk/ols/ontologies/envo
    | Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple terms can be separated by pipes, e.g.: shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]
* Field Format: free text
* Expected value: EnvO
* Value syntax: {term}
* Example: shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]

**environment (material)** : ``env_material``

* Definition:
    | The environmental material level refers to the material that was displaced by the sample, or material in which a sample was embedded, prior to the sampling event. Environmental matter terms are generally mass nouns. Examples include: air, soil, or water. EnvO (v1.53) terms listed under environmental matter can be found from the link: http://www.environmentontology.org/Browse-EnvO or https://www.ebi.ac.uk/ols/ontologies/envo
    | Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environmental material [ENVO:00010483]. Multiple terms can be separated by pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]
* Field Format: free text
* Expected value: EnvO
* Value syntax: {term}
* Example: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]

**geographic location** : ``geo_loc_name``

* Definition: The geographical origin of the sample as defined by the country or sea name followed by specific region name. Country or sea names should be chosen from the INSDC country list (http://insdc.org/country.html), or the GAZ ontology (v 1.512) (http://purl.bioontology.org/ontology/GAZ). Use a colon to separate the country or ocean from more detailed information about the location, eg "Canada: Vancouver" or "Germany: halfway down Zugspitze, Alps"
* Field Format: restricted text
* Expected value: country or sea name (INSDC or GAZ):region(GAZ):specific location name
* Value syntax: {term}:{term}:{text}
* Example: Germany\:Sylt\:Hausstrand

**latitude and longitude** : ``lat_lon``

* Definition: The geographical coordinates of the location where the sample was collected. The values should be reported in decimal degrees and in WGS84 system. Specify as degrees latitude and longitude in format "d[d.dddd] N\|S d[dd.dddd] W\|E", eg, 38.98 N 77.11 W
* Field Format: restricted text
* Expected value: decimal degrees
* Value syntax: {float} {float}
* Example: 38.98 N 77.11 W

organism characteristics
~~~~~~~~~~~~~~~~~~~~~~~~

**age** : ``age``

* Definition: Age at the time of sampling; relevant scale depends on species and study, e.g. could be seconds for amoebae or centuries for trees.
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: centuries,days,decades,hours,minutes,months,seconds,weeks,years

**beta-lactamase family** : ``beta_lactamase_family``

* Definition: Specify the beta-lactamase family for this gene.
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'ACC', 'ACT', 'ADC', 'BEL', 'CARB', 'CBP', 'CFE', 'CMY', 'CTX-M', 'DHA', 'FOX', 'GES', 'GIM', 'KPC', 'IMI', 'IMP', 'IND', 'LAT', 'MIR', 'MOX', 'NDM', 'OXA', 'PER', 'PDC', 'SHV', 'SME', 'TEM', 'VEB', 'VIM', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**carbapenemase** : ``carbapenemase``

* Definition: Does the enzyme exhibit carbapenemase activity? If the enzyme does exhibit carbapenemase activity, the response should be "yes", otherwise "no".
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'yes', 'no', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**ecotype** : ``ecotype``

* Definition: A population within a given species displaying genetically based, phenotypic traits that reflect adaptation to a local habitat, e.g., Columbia

**EDTA inhibitor tested** : ``edta_inhibitor_tested``

* Definition: Was carbapenemase activity tested in the presence of EDTA? If carbapenemase activity was tested in the presence of EDTA, the response should be "yes", otherwise "no".
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'yes', 'no', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**encoded traits** : ``encoded_traits``

* Definition: traits like antibiotic resistance or xenobiotic degradation phenotypes for plasmids, converting genes for phage
* Field Format: free text
* Expected value: for plasmid: antibiotic resistance; for phage: converting genes
* Value syntax: {text}

**genotype** : ``genotype``

* Definition: Observed genotype.
* Field Format: free text
* Expected value: genotype
* Value syntax: {text}

**height or length** : ``height_or_length``

* Definition: Measurement of height or length.
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}

**karyotype** : ``karyotype``

**mating type** : ``mating_type``

**observed biotic relationship** : ``biotic_relationship``

* Definition: Is it free-living or in a host and if the latter what type of relationship is observed
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'free living', 'parasite', 'commensal', 'symbiont', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**organism** : ``organism``

* Definition: The most descriptive organism name for this sample (to the species, if relevant).

**pathogenicity** : ``pathogenicity``

* Definition: To what is the entity pathogenic, for instance plant, fungi, bacteria
* Field Format: free text
* Expected value: CV
* Value syntax: {term}
* Example: human, animal, plant, fungi, bacteria

**pathotype** : ``pathotype``

* Definition: Some bacterial specific pathotypes (example Eschericia coli - STEC, UPEC).

**phenotype** : ``phenotype``

* Definition: Phenotype of sampled organism. For Phenotypic Quality Ontology (PATO) (v1.269) terms, please see http://purl.bioontology.org/ontology/PATO
* Field Format: free text
* Expected value: PATO
* Value syntax: {term}

**population** : ``population``

* Definition: for human: a collection of humans; for plants: filial generation, number of progeny, genetic structure

**race** : ``race``

**relationship to oxygen** : ``rel_to_oxygen``

* Definition: Is this organism an aerobe, anaerobe? Please note that aerobic and anaerobic are valid descriptors for microbial environments
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'aerobe', 'anaerobe', 'facultative', 'microaerophilic', 'microanaerobe', 'obligate aerobe', 'obligate anaerobe', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**serotype** : ``serotype``

* Definition: Taxonomy below subspecies; a variety (in bacteria, fungi or virus) usually based on its antigenic properties. Same as serovar and serogroup. e.g. serotype="H1N1" in Influenza A virus CY098518.

**serovar** : ``serovar``

* Definition: Taxonomy below subspecies; a variety (in bacteria, fungi or virus) usually based on its antigenic properties. Same as serotype and serogroup. Sometimes used as species identifier in bacteria with shaky taxonomy, e.g. Leptospira, serovar saopaolo S76607 (65357 in Entrez).

**sex** : ``sex``

* Definition: Physical sex of sampled organism.
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'male', 'female', 'pooled male and female', 'neuter', 'hermaphrodite', 'intersex', 'not determined', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**strain** : ``strain``

* Definition: Microbial or eukaryotic strain name.

**subgroup** : ``subgroup``

* Definition: Taxonomy below subspecies; sometimes used in viruses to denote subgroups taken from a single isolate.

**subspecific genetic lineage** : ``subspecf_gen_lin``

* Definition: This should provide further information about the genetic distinctness of this lineage by recording additional information i.e biovar, serovar, serotype, biovar, or any relevant genetic typing schemes like Group I plasmid. It can also contain alternative taxonomic information
* Field Format: free text
* Expected value: genetic lineage below lowest rank of NCBI taxonomy, which is subspecies, e.g. serovar, biotype, ecotype
* Value syntax: {text}

**subtype** : ``subtype``

* Definition: Used as classifier in viruses (e.g. HIV type 1, Group M, Subtype A).

**taxonomy ID** : ``tax_id``

* Definition: The Taxonomy ID indicates the taxonomic classification of the sample (e.g. 9606 for human). For metagenomic samples please browse the http://www.ebi.ac.uk/ena/data/view/Taxon:408169 for a suitable taxon. For a previously unsequenced organisms please contact datasubs@cngb.org for the provision of a new Taxonomy ID.

**trophic level** : ``trophic_level``

* Definition: Trophic levels are the feeding position in a food chain. Microbes can be a range of producers (e.g. chemolithotroph)
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'autotroph', 'carboxydotroph', 'chemoautotroph', 'chemoheterotroph', 'chemolithoautotroph', 'chemolithotroph', 'chemoorganoheterotroph', 'chemoorganotroph', 'chemosynthetic', 'chemotroph', 'copiotroph', 'diazotroph', 'facultative', 'heterotroph', 'lithoautotroph', 'lithoheterotroph', 'lithotroph', 'methanotroph', 'methylotroph', 'mixotroph', 'obligate', 'chemoautolithotroph', 'oligotroph', 'organoheterotroph', 'organotroph', 'photoautotroph', 'photoheterotroph', 'photolithoautotroph', 'photolithotroph', 'photosynthetic', 'phototroph', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

sample collection
~~~~~~~~~~~~~~~~~

**biomass** : ``biomass``

* Definition: amount of biomass; should include the name for the part of biomass measured, e.g. microbial, total. can include multiple measurements
* Field Format: free text
* Expected value: biomass type;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: ton, kilogram, gram

**biomaterial provider** : ``biomaterial_provider``

* Definition: Name and address of the lab or PI, or a culture collection identifier.
* Field Format: free text

**birth date** : ``birth_date``

**birth location** : ``birth_location``

**breed** : ``breed``

* Definition: breed name - chiefly used in domesticated animals or plants.

**breeding history** : ``breeding_history``

**breeding method** : ``breeding_method``

**cell line** : ``cell_line``

* Definition: Name of the cell line.

**cell subtype** : ``cell_subtype``

**cell type** : ``cell_type``

* Definition: Type of cell of the sample or from which the sample was obtained.

**collected by** : ``collected_by``

* Definition: Name of persons or institute who collected the sample.
* Field Format: free text

**cultivar** : ``cultivar``

* Definition: cultivar name - cultivated variety of plant.

**culture collection** : ``culture_collection``

* Definition: Name of source institute and unique culture identifier. See the description for the proper format and list of allowed institutes, http://www.insdc.org/controlled-vocabulary-culturecollection-qualifier
* Field Format: restricted text
* Expected value: <institution-code>:[<collection-code>:]<culture_id>
* Value syntax: {term}:{term}:{text}
* Example: ATCC:26370

**death date** : ``death_date``

**density** : ``density``

* Definition: density of sample
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: gram per cubic meter

**development stage** : ``dev_stage``

* Definition: Developmental stage at the time of sampling.

**disease** : ``disease``

* Definition: List of diseases diagnosed; can include multiple diagnoses. the value of the field depends on host; for humans the terms should be chosen from DO (Disease Ontology), free text for non-human. For DO terms, please see https://www.ebi.ac.uk/ols/ontologies/symp
* Field Format: free text
* Expected value: disease name or DO
* Value syntax: {term}

**disease stage** : ``disease_stage``

* Definition: Stage of disease at the time of sampling.

**growth protocol** : ``growth_protocol``

**host body product** : ``host_body_product``

* Definition: substance produced by the host, e.g. stool, mucus, where the sample was obtained from. For Foundational Model of Anatomy Ontology (FMA) (v 3.1) terms, please see http://purl.bioontology.org/ontology/FMA
* Field Format: free text
* Expected value: FMA
* Value syntax: {term}

**host dry mass** : ``host_dry_mass``

* Definition: measurement of dry mass
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: kilogram, gram

**host wet mass** : ``host_wet_mass``

* Definition: measurement of wet mass
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: kilogram, gram

**identified by** : ``identified_by``

* Definition: Name of the taxonomist who identified the specimen.

**isolate** : ``isolate``

* Definition: Identification or description of the specific individual from which this sample was obtained.

**isolation and growth condition** : ``isol_growth_condt``

* Definition: Publication reference in the form of pubmed ID (pmid), digital object identifier (doi) or url for isolation and growth condition specifications of the organism/material.
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**isolation source** : ``isolation_source``

* Definition: Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived.

**medical history performed** : ``medic_hist_perform``

* Definition: whether full medical history was collected
* Field Format: text choice
* Expected value: medical history status
* Value syntax: {boolean}

**microbial biomass** : ``microbial_biomass``

* Definition: the part of the organic matter in the soil that constitutes living microorganisms smaller than 5-10 micrometer. If you keep this, you would need to have correction factors used for conversion to the final units, which should be mg C (or N)/kg soil).
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: ton, kilogram, gram per kilogram soil

**microbial biomass method** : ``microbial_biomass_meth``

* Definition: reference or method used in determining microbial biomass
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**Omics Observatory ID** : ``omics_observ_id``

* Definition: A unique identifier of the omics-enabled observatory (or comparable time series) your data derives from. This identifier should be provided by the OMICON ontology; if you require a new identifier for your time series, contact the ontology's developers. Information is available here: https://github.com/GLOMICON/omicon. This field is only applicable to records which derive from an omics time-series or observatory.
* Field Format: free text
* Expected value: OMICON
* Value syntax: {term}

**organism count** : ``organism_count``

* Definition: total count of any organism per gram or volume of sample, should include name of organism followed by count; can include multiple organism counts
* Field Format: free text
* Expected value: organism name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: number of organism per cubic meter

**oxygenation status of sample** : ``oxy_stat_samp``

* Definition: oxygenation status of sample
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'aerobic', 'anaerobic', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**passage history** : ``passage_history``

* Definition: Number of passages and passage method.

**plant product** : ``plant_product``

* Definition: substance produced by the plant, where the sample was obtained from
* Field Format: free text
* Expected value: product name
* Value syntax: {text}

**propagation** : ``propagation``

* Definition: This field is specific to different taxa. For phages: lytic/lysogenic/temperate/obligately lytic, for plasmids: incompatibility group, for eukaryote: asexual/sexual
* Field Format: free text
* Expected value: for virus: lytic, lysogenic, temperate, obligately lytic; for plasmid: incompatibility group; for eukaryote: asexual, sexual) [CV]
* Value syntax: {term}

**sample collection device or method** : ``samp_collect_device``

* Definition: The method or device employed for collecting the sample
* Field Format: free text
* Expected value: type name
* Value syntax: {text}
* Example: biopsy, niskin bottle, push core

**sample material processing** : ``samp_mat_process``

* Definition: Any processing applied to the sample during or after retrieving the sample from environment. This field accepts OBI, for a browser of OBI (v 2013-10-25) terms please see http://purl.bioontology.org/ontology/OBI
* Field Format: free text
* Expected value: text or OBI
* Value syntax: {text|term}
* Example: filtering of seawater, storing samples in ethanol

**sample name** : ``sample_name``

* Definition: Sample name is a name that you choose for the sample. It can have any format, but we suggest that you make it concise, unique and consistent within your lab, and as informative as possible. Every Sample Name from a single Submitter must be unique.

**sample size** : ``samp_size``

* Definition: Amount or size of sample (volume, mass or area) that was collected
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: liter,gram,kilogram,square meter,cubic meter

**sample size sorting method** : ``samp_sort_meth``

* Definition: method by which samples are sorted
* Field Format: free text
* Expected value: description of method
* Value syntax: {text}
* Example: open face filter collecting total suspended particles, "prefilter to remove particles larger than X micrometers in diameter, where common values of X would be 10 and 2.5 full size sorting in a cascade impactor"

**sample storage duration** : ``samp_store_dur``

* Definition: duration for which sample was stored
* Field Format: restricted text
* Expected value: time interval
* Value syntax: {interval}

**sample storage location** : ``samp_store_loc``

* Definition: location at which sample was stored, usually name of a specific freezer/room
* Field Format: free text
* Expected value: location name
* Value syntax: {text}

**sample storage temperature** : ``samp_store_temp``

* Definition: temperature at which sample was stored, e.g. -80 degree Celsius
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree Celsius

**sample title** : ``sample_title``

* Definition: The sample title is a short, preferably a single sentence, description of the sample.

**sample type** : ``sample_type``

* Definition: Sample type, such as cell culture, mixed culture, tissue sample, whole organism, single cell, metagenomic assembly.

**sieving** : ``sieving``

* Definition: collection design of pooled samples and/or sieve size and amount of sample sieved
* Field Format: free text
* Expected value: design name and/or size;amount
* Value syntax: {{text}|{float} {unit}};{float} {unit}

**size fraction selected** : ``size_frac``

* Definition: Filtering pore size used in sample preparation, e.g., 0-0.22 micrometer
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float}-{float} {unit}

**source material identifiers** : ``source_material_id``

* Definition: A unique identifier assigned to a material sample (as defined by http://rs.tdwg.org/dwc/terms/materialSampleID, and as opposed to a particular digital record of a material sample) used for extracting nucleic acids, and subsequent sequencing. The identifier can refer either to the original material collected or to any derived sub-samples.
* Field Format: free text
* Expected value: for cultures of microorganisms: identifiers for two culture collections; for other material a unique arbitrary identifer
* Value syntax: {text}
* Example: xatc123, ark:/2154/R2

**source of UViGs** : ``source_uvig``

* Definition: Type of dataset from which the UViG was obtained
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'metagenome (not viral targeted)', 'viral fraction metagenome (virome)', 'sequence-targeted metagenome', 'metatranscriptome (not viral targeted)', 'viral fraction RNA metagenome (RNA virome)', 'sequence-targeted RNA metagenome', 'microbial single amplified genome (SAG)', 'viral single amplified genome (vSAG)', 'isolate microbial genome', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**specimen voucher** : ``specimen_voucher``

* Definition: Identifier for the physical specimen. Use format: "[<institution-code>:[<collection-code>:]]<specimen_id>", eg, "UAM:Mamm:52179". Intended as a reference to the physical specimen that remains after it was analyzed. If the specimen was destroyed in the process of analysis, electronic images (e-vouchers) are an adequate substitute for a physical voucher specimen. Ideally the specimens will be deposited in a curated museum, herbarium, or frozen tissue collection, but often they will remain in a personal or laboratory collection for some time before they are deposited in a curated collection. There are three forms of specimen_voucher qualifiers. If the text of the qualifier includes one or more colons it is a 'structured voucher'. Structured vouchers include institution-codes (and optional collection-codes) taken from a controlled vocabulary maintained by the INSDC that denotes the museum or herbarium collection where the specimen resides, please visit: http://www.insdc.org/controlled-vocabulary-specimenvoucher-qualifier.
* Field Format: restricted text
* Expected value: [<institution-code>:[<collection-code>:]]<specimen_id>
* Value syntax: {term}:{term}:{text}
* Example: UAM\:Mamm\:52179

**storage conditions** : ``store_cond``

* Definition: explain how and for how long the soil sample was stored before DNA extraction.
* Field Format: free text
* Expected value: storage condition type;duration
* Value syntax: {text};{period}

**stud book number** : ``stud_book_number``

**tissue** : ``tissue``

* Definition: Type of tissue the sample was taken from.

**urine collection method** : ``urine_collect_meth``

* Definition: specification of urine collection method
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'clean catch', 'catheter', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**virus enrichment approach** : ``virus_enrich_appr``

* Definition: Approach used to enrich the sample for viruses, if any. If more than one approach was used, include multiple ‘virus_enrich_appr’ fields.
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'filtration', 'ultrafiltration', 'centrifugation', 'ultracentrifugation', 'PEG Precipitation', 'FeCl Precipitation', 'CsCl density gradient', 'DNAse', 'RNAse', 'targeted sequence capture', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']
* Example: filtration + FeCl Precipitation + ultracentrifugation + DNAse

host description
~~~~~~~~~~~~~~~~

**drug usage** : ``drug_usage``

* Definition: any drug used by subject and the frequency of usage; can include multiple drugs used
* Field Format: free text
* Expected value: drug name;frequency
* Value syntax: {text};{integer}/[year\|month\|week\|day\|hour]

**ethnicity** : ``ethnicity``

* Definition: Ethnicity of the subject.
* Field Format: free text
* Expected value: IHMC code or free text
* Value syntax: {integer|text}

**health state ** : ``health_state``

* Definition: Health or disease status of specific host at time of collection. This field accepts PATO (v1.269) terms, for a browser please see http://purl.bioontology.org/ontology/PATO
* Field Format: free text
* Expected value: PATO
* Value syntax: {term}

**host** : ``host``

* Definition: If there is a host involved, please provide its taxid (or environmental if not actually isolated from the dead or alive host - i.e. pathogen could be isolated from a swipe of a bench etc) and report whether it is a laboratory or natural host).
* Field Format: free text
* Expected value: host taxid, unknown, environmental
* Example: Homo sapiens

**host age** : ``host_age``

* Definition: age of host at the time of sampling; relevant scale depends on species and study, e.g. could be seconds for amoebae or centuries for trees
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: centuries,days,decades,hours,minutes,months,seconds,weeks,years

**host body habitat** : ``host_body_habitat``

* Definition: original body habitat where the sample was obtained from
* Field Format: free text
* Expected value: FMA
* Value syntax: {term}

**host body mass index** : ``host_body_mass_index``

* Definition: body mass index of the host, calculated as weight/(height)squared
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: kilogram per square meter

**host body temperature** : ``host_body_temp``

* Definition: core body temperature of the host when sample was collected
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree Celsius

**host color** : ``host_color``

* Definition: the color of host
* Field Format: free text
* Expected value: color
* Value syntax: {text}

**host description** : ``host_description``

* Definition: Additional information not included in other defined vocabulary fields.

**host disease outcome** : ``host_disease_outcome``

* Definition: Final outcome of disease, e.g., death, chronic disease, recovery.

**host disease stage** : ``host_disease_stage``

* Definition: Stage of disease at the time of sampling.

**host growth conditions** : ``host_growth_cond``

* Definition: literature reference giving growth conditions of the host
* Field Format: free text
* Expected value: PMID,DOI,url or free text
* Value syntax: {PMID\|DOI\|URL}

**host health state** : ``host_health_state``

* Definition: Information regarding health state of the individual sampled at the time of sampling.
* Field Format: free text

**host height** : ``host_height``

* Definition: the height of subject
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: centimeter, millimeter, meter

**host infra-specific name** : ``host_infra_specific_name``

* Definition: taxonomic information about the host below subspecies level
* Field Format: free text
* Expected value: name
* Value syntax: {text}

**host infra-specific rank** : ``host_infra_specific_rank``

* Definition: taxonomic rank information about the host below subspecies level, such as variety, form, rank etc.
* Field Format: free text
* Expected value: rank
* Value syntax: {text}

**host length** : ``host_length``

* Definition: the length of subject
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: centimeter, millimeter, meter

**host life stage** : ``host_life_stage``

* Definition: description of life stage of host
* Field Format: free text
* Expected value: stage
* Value syntax: {text}

**host occupation** : ``host_occupation``

* Definition: most frequent job performed by subject
* Field Format: text choice
* Expected value: IHMC code
* Value syntax: {integer}

**host phenotype** : ``host_phenotype``

* Definition: phenotype of host. For Phenotypic Quality Ontology (PATO) (v1.269) terms, please see http://purl.bioontology.org/ontology/PATO
* Field Format: free text
* Expected value: PATO
* Value syntax: {term}

**host sex** : ``host_sex``

* Definition: physical sex of the host
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'male', 'female', 'pooled male and female', 'neuter', 'hermaphrodite', 'intersex', 'not determined', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**host shape** : ``host_shape``

* Definition: morphological shape of host
* Field Format: free text
* Expected value: shape
* Value syntax: {text}

**host subject id** : ``host_subject_id``

* Definition: a unique identifier by which each subject can be referred to, de-identified, e.g. #131
* Field Format: free text
* Expected value: unique identifier
* Value syntax: {text}

**host substrate** : ``host_substrate``

* Definition: the growth substrate of the host
* Field Format: free text
* Expected value: substrate name
* Value syntax: {text}

**host taxonomy ID** : ``host_taxid``

* Definition: NCBI taxonomy ID of the host, e.g. 9606
* Field Format: restricted text
* Expected value: NCBI taxon identifier
* Value syntax: {integer}

**host tissue sampled** : ``host_tissue_sampled``

* Definition: name of body site where the sample was obtained from, such as a specific organ or tissue (tongue, lung etc...). For Foundational Model of Anatomy Ontology (FMA) (v 3.1) terms, please see http://purl.bioontology.org/ontology/FMA
* Field Format: free text
* Expected value: FMA
* Value syntax: {term}

**host total mass** : ``host_tot_mass``

* Definition: total mass of the host at collection, the unit depends on host
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: kilogram, gram

**IHMC medication code** : ``ihmc_medication_code``

* Definition: can include multiple medication codes
* Field Format: text choice
* Expected value: IHMC code
* Value syntax: {integer}

**lab host** : ``lab_host``

* Definition: Scientific name and description of the laboratory host used to propagate the source organism or material from which the sample was obtained, e.g., Escherichia coli DH5a, or Homo sapiens HeLa cells.

**plant body site** : ``plant_body_site``

* Definition: name of body site that the sample was obtained from. For Plant Ontology (PO) (v 20) terms, see http://purl.bioontology.org/ontology/PO
* Field Format: free text
* Expected value: PO
* Value syntax: {term}

host details
~~~~~~~~~~~~

**amniotic fluid color** : ``amniotic_fluid_color``

* Definition: specification of the color of the amniotic fluid sample
* Field Format: free text
* Expected value: color
* Value syntax: {text}

**birth control** : ``birth_control``

* Definition: specification of birth control medication used
* Field Format: free text
* Expected value: medication name
* Value syntax: {text}

**dominant hand** : ``dominant_hand``

* Definition: dominant hand of the subject
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'left', 'right', 'ambidextrous', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**douche** : ``douche``

* Definition: date of most recent douche
* Field Format: free text
* Expected value: timestamp
* Value syntax: {timestamp}

**foetal health status** : ``foetal_health_stat``

* Definition: specification of foetal health status, should also include abortion
* Field Format: free text
* Expected value: health status
* Value syntax: {text}

**gestation state** : ``gestation_state``

* Definition: specification of the gestation state
* Field Format: free text
* Expected value: gestation state
* Value syntax: {text}

**gravidity** : ``gravidity``

* Definition: information about treatment involving use of gravity factor to study various types of responses in presence, absence or modified levels of gravity; can include multiple treatments
* Field Format: free text
* Expected value: gravity factor value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: meter per square second, g

**host blood pressure diastolic** : ``host_blood_press_diast``

* Definition: resting diastolic blood pressure of the host, measured as mm mercury
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millimeter mercury

**host blood pressure systolic** : ``host_blood_press_syst``

* Definition: resting systolic blood pressure of the host, measured as mm mercury
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millimeter mercury

**host diet** : ``host_diet``

* Definition: type of diet depending on the host, for animals omnivore, herbivore etc., for humans high-fat, meditteranean etc.; can include multiple diet types
* Field Format: free text
* Expected value: diet type
* Value syntax: {text}

**host family relationship** : ``host_family_relationship``

* Definition: relationships to other hosts in the same study; can include multiple relationships
* Field Format: free text
* Expected value: relationship type;arbitrary identifier
* Value syntax: {text};{text}

**host genotype** : ``host_genotype``

* Definition: observed genotype
* Field Format: free text
* Expected value: genotype
* Value syntax: {text}

**host last meal** : ``host_last_meal``

* Definition: content of last meal and time since feeding; can include multiple values
* Field Format: free text
* Expected value: content;time interval
* Value syntax: {text};{period}

**host pulse** : ``host_pulse``

* Definition: resting pulse of the host, measured as beats per minute
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: beats per minute

**HRT** : ``hrt``

* Definition: whether subject had hormone replacement theraphy, and if yes start date
* Field Format: free text
* Expected value: timestamp
* Value syntax: {timestamp}

**hysterectomy** : ``hysterectomy``

* Definition: specification of whether hysterectomy was performed
* Field Format: text choice
* Expected value: hysterectomy status
* Value syntax: {boolean}
* Example: yes,no

**major diet change in last six months** : ``diet_last_six_month``

* Definition: specification of major diet changes in the last six months, if yes the change should be specified
* Field Format: free text
* Expected value: diet change;current diet
* Value syntax: {boolean};{text}

**maternal health status** : ``maternal_health_stat``

* Definition: specification of the maternal health status
* Field Format: free text
* Expected value: health status
* Value syntax: {text}

**menarche** : ``menarche``

* Definition: date of most recent menstruation
* Field Format: free text
* Expected value: timestamp
* Value syntax: {timestamp}

**menopause** : ``menopause``

* Definition: date of onset of menopause
* Field Format: free text
* Expected value: timestamp
* Value syntax: {timestamp}

**pregnancy** : ``pregnancy``

* Definition: date due of pregnancy
* Field Format: free text
* Expected value: timestamp
* Value syntax: {timestamp}

**sexual activity** : ``sexual_act``

* Definition: current sexual partner and frequency of sex
* Field Format: free text
* Expected value: partner sex;frequency
* Value syntax: {text}

**smoker** : ``smoker``

* Definition: specification of smoking status
* Field Format: text choice
* Expected value: smoking status
* Value syntax: ['', 'ex-smoker', 'non-smoker', 'smoker', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**special diet** : ``special_diet``

* Definition: specification of special diet; can include multiple special diets
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'low carb', 'reduced calorie', 'vegetarian', 'other(to be specified)', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**time since last toothbrushing** : ``time_last_toothbrush``

* Definition: specification of the time since last toothbrushing
* Field Format: restricted text
* Expected value: timestamp
* Value syntax: {timestamp}

**time since last wash** : ``time_since_last_wash``

* Definition: specification of the time since last wash
* Field Format: restricted text
* Expected value: timestamp
* Value syntax: {timestamp}

**travel outside the country in last six months** : ``travel_out_six_month``

* Definition: specification of the countries travelled in the last six months; can include multiple travels
* Field Format: free text
* Expected value: country name
* Value syntax: {text}

**twin sibling presence** : ``twin_sibling``

* Definition: specification of twin sibling presence
* Field Format: text choice
* Expected value: presence status
* Value syntax: {boolean}

**weight loss in last three months** : ``weight_loss_3_month``

* Definition: specification of weight loss in the last three months, if yes should be further specified to include amount of weight loss
* Field Format: free text
* Expected value: weight loss specification;measurement value
* Value syntax: {boolean};{float} {unit}
* Preferred unit: kilogram, gram

host disorder
~~~~~~~~~~~~~

**blood disorder** : ``blood_blood_disord``

* Definition: history of blood disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**dermatology disorder** : ``dermatology_disord``

* Definition: history of dermatology disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**gastrointestinal tract disorder** : ``gastrointest_disord``

* Definition: history of gastrointestinal tract disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**gynecological disorder** : ``gynecologic_disord``

* Definition: history of gynecological disorders; can include multiple disorders
* Field Format: free text
* Expected value: gynecological disorder
* Value syntax: {text}

**host disease** : ``host_disease``

* Definition: list of diseases with which the host has been diagnosed; can include multiple diagnoses. the value of the field depends on host; for humans the terms should be chosen from DO (Disease Ontology) at http://www.disease-ontology.org, other hosts are free text
* Field Format: free text
* Expected value: disease name or DO
* Value syntax: {term}

**host HIV status** : ``host_hiv_stat``

* Definition: HIV status of subject, if yes HAART initiation status should also be indicated as [YES or NO]
* Field Format: text choice
* Expected value: HIV status;HAART initiation status
* Value syntax: {boolean};{boolean}

**kidney disorder** : ``kidney_disord``

* Definition: history of kidney disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**liver disorder** : ``liver_disord``

* Definition: history of liver disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**nose/mouth/teeth/throat disorder** : ``nose_mouth_teeth_throat_disord``

* Definition: history of nose/mouth/teeth/throat disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**nose-throat disorder** : ``nose_throat_disord``

* Definition: history of nose-throat disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**pulmonary disorder** : ``pulmonary_disord``

* Definition: history of pulmonary disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**urogenital disorder** : ``urogenit_disord``

* Definition: history of urogenital disorders, can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

**urogenital tract disorder** : ``urogenit_tract_disor``

* Definition: history of urogenitaltract disorders; can include multiple disorders
* Field Format: free text
* Expected value: disorder name
* Value syntax: {text}

bioreactor
~~~~~~~~~~

**biochemical oxygen demand** : ``biochem_oxygen_dem``

* Definition: a measure of the relative oxygen-depletion effect of a waste contaminant
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter

**chemical oxygen demand** : ``chem_oxygen_dem``

* Definition: a measure of the relative oxygen-depletion effect of a waste contaminant
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter

**pre-treatment** : ``pre_treatment``

* Definition: the process of pre-treatment removes materials that can be easily collected from the raw wastewater
* Field Format: free text
* Expected value: pre-treatment type
* Value syntax: {text}

**primary treatment** : ``primary_treatment``

* Definition: the process to produce both a generally homogeneous liquid capable of being treated biologically and a sludge that can be separately treated or processed
* Field Format: free text
* Expected value: primary treatment type
* Value syntax: {text}

**reactor type** : ``reactor_type``

* Definition: anaerobic digesters can be designed and engineered to operate using a number of different process configurations, as batch or continuous, mesophilic, high solid or low solid, and single stage or multistage
* Field Format: free text
* Expected value: reactor type name
* Value syntax: {text}

**secondary treatment** : ``secondary_treatment``

* Definition: the process for substantially degrading the biological content of the sewage
* Field Format: free text
* Expected value: secondary treatment type
* Value syntax: {text}

**sludge retention time** : ``sludge_retent_time``

* Definition: the time activated sludge remains in reactor
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: days,hours,minutes,weeks

**tertiary treatment** : ``tertiary_treatment``

* Definition: the process providing a final treatment stage to raise the effluent quality before it is discharged to the receiving environment
* Field Format: free text
* Expected value: tertiary treatment type
* Value syntax: {text}

**treatment** : ``treatment``

concentration measurement
~~~~~~~~~~~~~~~~~~~~~~~~~

**alkyl diethers** : ``alkyl_diethers``

* Definition: concentration of alkyl diethers
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: mole per liter

**aminopeptidase activity** : ``aminopept_act``

* Definition: measurement of aminopeptidase activity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: mole per liter per hour

**ammonium** : ``ammonium``

* Definition: concentration of ammonium
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**bacterial carbon production** : ``bacteria_carb_prod``

* Definition: measurement of bacterial carbon production
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: nanogram per hour

**bacterial production** : ``bac_prod``

* Definition: bacterial production in the water column measured by isotope uptake
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per cubic meter per day

**bacterial respiration** : ``bac_resp``

* Definition: measurement of bacterial respiration in the water column
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per cubic meter per day

**bishomohopanol** : ``bishomohopanol``

* Definition: concentration of bishomohopanol
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter, microgram per gram

**bromide** : ``bromide``

* Definition: concentration of bromide
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: parts per million

**calcium** : ``calcium``

* Definition: concentration of calcium
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter, micromole per liter, parts per million

**carbon dioxide** : ``carb_dioxide``

* Definition: carbon dioxide (gas) amount or concentration at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**carbon monoxide** : ``carb_monoxide``

* Definition: carbon monoxide (gas) amount or concentration at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**carbon/nitrogen ratio** : ``carb_nitro_ratio``

* Definition: ratio of amount or concentrations of carbon to nitrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}

**chloride** : ``chloride``

* Definition: concentration of chloride
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter

**chlorophyll** : ``chlorophyll``

* Definition: concentration of chlorophyll
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per cubic meter, microgram per liter

**diether lipids** : ``diether_lipids``

* Definition: concentration of diether lipids; can include multiple types of diether lipids
* Field Format: restricted text
* Expected value: diether lipid name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: nanogram per liter

**dissolved carbon dioxide** : ``diss_carb_dioxide``

* Definition: concentration of dissolved carbon dioxide
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**dissolved hydrogen** : ``diss_hydrogen``

* Definition: concentration of dissolved hydrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**dissolved inorganic carbon** : ``diss_inorg_carb``

* Definition: dissolved inorganic carbon concentration
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**dissolved inorganic nitrogen** : ````

* Definition: concentration of dissolved inorganic nitrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**dissolved inorganic phosphorus** : ``diss_inorg_phosp``

* Definition: concentration of dissolved inorganic phosphorus
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**dissolved organic carbon** : ``diss_org_carb``

* Definition: concentration of dissolved organic carbon
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**dissolved organic nitrogen** : ``diss_org_nitro``

* Definition: dissolved organic nitrogen concentration measured as; total dissolved nitrogen - NH4 - NO3 - NO2
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter, milligram per liter

**dissolved oxygen** : ``diss_oxygen``

* Definition: concentration of dissolved oxygen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per kilogram

**efficiency percent** : ``efficiency_percent``

* Definition: percentage of volatile solids removed from the anaerobic digestor
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage

**emulsions** : ``emulsions``

* Definition: amount or concentration of substances such as paints, adhesives, mayonnaise, hair colorants, emulsified oils, etc.; can include multiple emulsion types
* Field Format: restricted text
* Expected value: emulsion name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram per liter

**gaseous substances** : ``gaseous_substances``

* Definition: amount or concentration of substances such as hydrogen sulfide, carbon dioxide, methane, etc.; can include multiple substances
* Field Format: restricted text
* Expected value: gaseous substance name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: micromole per liter

**glucosidase activity** : ``glucosidase_act``

* Definition: measurement of glucosidase activity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: mol per liter per hour

**inorganic particles** : ``inorg_particles``

* Definition: concentration of particles such as sand, grit, metal particles, ceramics, etc.; can include multiple particles
* Field Format: restricted text
* Expected value: inorganic particle name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: mole per liter, milligram per liter

**magnesium** : ``magnesium``

* Definition: concentration of magnesium
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: mole per liter, milligram per liter, parts per million

**methane** : ``methane``

* Definition: methane (gas) amount or concentration at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**n-alkanes** : ``n_alkanes``

* Definition: concentration of n-alkanes; can include multiple n-alkanes
* Field Format: restricted text
* Expected value: n-alkane name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: micromole per liter

**nitrate** : ``nitrate``

* Definition: concentration of nitrate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**nitrite** : ``nitrite``

* Definition: concentration of nitrite
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**nitrogen** : ``nitro``

* Definition: concentration of nitrogen (total)
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**organic carbon** : ``org_carb``

* Definition: concentration of organic carbon
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**organic matter** : ``org_matter``

* Definition: concentration of organic matter
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**organic nitrogen** : ``org_nitro``

* Definition: concentration of organic nitrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**organic particles** : ``org_particles``

* Definition: concentration of particles such as faeces, hairs, food, vomit, paper fibers, plant material, humus, etc.
* Field Format: restricted text
* Expected value: particle name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram per liter

**oxygen** : ``oxygen``

* Definition: oxygen (gas) amount or concentration at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter, parts per million

**particulate organic carbon** : ``part_org_carb``

* Definition: concentration of particulate organic carbon
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**particulate organic nitrogen** : ``part_org_nitro``

* Definition: concentration of particulate organic nitrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**petroleum hydrocarbon** : ``petroleum_hydrocarb``

* Definition: concentration of petroleum hydrocarbon
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**phaeopigments** : ``phaeopigments``

* Definition: concentration of phaeopigments; can include multiple phaeopigments
* Field Format: restricted text
* Expected value: phaeopigment name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: milligram per cubic meter

**phosphate** : ``phosphate``

* Definition: concentration of phosphate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**phospholipid fatty acid** : ``phosplipid_fatt_acid``

* Definition: concentration of phospholipid fatty acids; can include multiple values
* Field Format: restricted text
* Expected value: phospholipid fatty acid name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: mole per gram, mole per liter

**potassium** : ``potassium``

* Definition: concentration of potassium
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: parts per million

**primary production** : ``primary_prod``

* Definition: measurement of primary production
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per cubic meter per day, gram per square meter per day

**redox potential** : ``redox_potential``

* Definition: redox potential, measured relative to a hydrogen cell, indicating oxidation or reduction potential
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millivolt

**respirable particulate matter** : ``resp_part_matter``

* Definition: concentration of substances that remain suspended in the air, and comprise mixtures of organic and inorganic substances (PM10 and PM2.5); can report multiple PM's by entering numeric values preceded by name of PM
* Field Format: restricted text
* Expected value: particulate matter name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: microgram per cubic meter

**salinity** : ``salinity``

* Definition: salinity measurement
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: practical salinity unit, percentage

**sample salinity** : ``samp_salinity``

* Definition: salinity of sample, i.e. measure of total salt concentration
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: practical salinity unit, percentage

**silicate** : ``silicate``

* Definition: concentration of silicate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**sodium** : ``sodium``

* Definition: sodium concentration
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: parts per million

**soluble inorganic material** : ``soluble_inorg_mat``

* Definition: concentration of substances such as ammonia, road-salt, sea-salt, cyanide, hydrogen sulfide, thiocyanates, thiosulfates, etc.
* Field Format: restricted text
* Expected value: soluble inorganic material name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram, microgram, mole per liter, gram per liter, parts per million

**soluble organic material** : ``soluble_org_mat``

* Definition: concentration of substances such as urea, fruit sugars, soluble proteins, drugs, pharmaceuticals, etc.
* Field Format: restricted text
* Expected value: soluble organic material name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram, microgram, mole per liter, gram per liter, parts per million

**soluble reactive phosphorus** : ``soluble_react_phosp``

* Definition: concentration of soluble reactive phosphorus
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**sulfate** : ``sulfate``

* Definition: concentration of sulfate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter, milligram per liter

**sulfide** : ``sulfide``

* Definition: concentration of sulfide
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter, milligram per liter

**suspended particulate matter** : ``suspend_part_matter``

* Definition: concentration of suspended particulate matter
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram per liter

**suspended solids** : ``suspend_solids``

* Definition: concentration of substances including a wide variety of material, such as silt, decaying plant and animal matter, etc,; can include multiple substances
* Field Format: restricted text
* Expected value: suspended solid name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram, microgram, mole per liter, gram per liter, part per million

**total carbon** : ``tot_carb``

* Definition: total carbon content
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**total dissolved nitrogen** : ``tot_diss_nitro``

* Definition: total dissolved nitrogen concentration, reported as nitrogen, measured by: total dissolved nitrogen = NH4 + NO3NO2 + dissolved organic nitrogen
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**total inorganic nitrogen** : ``tot_inorg_nitro``

* Definition: total inorganic nitrogen content
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter

**total N method** : ``tot_n_meth``

* Definition: reference or method used in determining the total nitrogen
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**total nitrogen** : ``tot_nitro``

* Definition: total nitrogen content of the sample
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter, micromole per liter

**total organic carbon** : ``tot_org_carb``

* Definition: Definition for soil: total organic C content of the soil units of g C/kg soil. Definition otherwise: total organic carbon content
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: gram Carbon per kilogram sample material

**total organic carbon method** : ``tot_org_c_meth``

* Definition: reference or method used in determining total organic carbon
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**total particulate carbon** : ``tot_part_carb``

* Definition: total particulate carbon content
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter, micromole per liter

**total phosphate** : ``tot_phosphate``

* Definition: total amount or concentration of phosphate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microgram per liter, micromole per liter

**total phosphorus** : ``tot_phosp``

* Definition: total phosphorus concentration, calculated by: total phosphorus = total dissolved phosphorus + particulate phosphorus. Can also be measured without filtering, reported as phosphorus
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per liter

**volatile organic compounds** : ``volatile_org_comp``

* Definition: concentration of carbon-based chemicals that easily evaporate at room temperature; can report multiple volatile organic compounds by entering numeric values preceded by name of compound
* Field Format: free text
* Expected value: volatile organic compound name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: microgram per cubic meter, parts per million

**water content** : ``water_content``

* Definition: water content measurement
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: gram per gram or cubic centimeter per cubic centimeter

**water content of soil** : ``water_content_soil``

* Definition: Water content (g/g or cm3/cm3).
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} [g/g|cm3/cm3]

**water content of soil method** : ``water_content_soil_meth``

* Definition: reference or method used in determining the water content of soil
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

geography
~~~~~~~~~

**profile position** : ``profile_position``

* Definition: cross-sectional position in the hillslope where sample was collected. sample area position in relation to surrounding areas
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'summit', 'shoulder', 'backslope', 'footslope', 'toeslope', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**slope aspect** : ``slope_aspect``

* Definition: the direction a slope faces. While looking down a slope use a compass to record the direction you are facing (direction or degrees); e.g., NW or 315 degrees. This measure provides an indication of sun and wind exposure that will influence soil temperature and evapotranspiration.
* Field Format: free text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree

**slope gradient** : ``slope_gradient``

* Definition: commonly called 'slope'. The angle between ground surface and a horizontal line (in percent). This is the direction that overland water would flow. This measure is usually taken with a hand level meter or clinometer
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage

internal environment
~~~~~~~~~~~~~~~~~~~~

**building occupancy type** : ``build_occup_type``

* Definition: the primary function for which a building or discrete part of a building is intended to be used
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'airport', 'agricultural', 'commercial', 'educational', 'government', 'health care', 'high rise', 'industrial', 'low rise', 'market', 'office', 'military', 'parking', 'residential', 'restaurant', 'school', 'sports complex', 'storage', 'religious', 'transport', 'wood framed', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**building setting** : ``building_setting``

* Definition: a location (geography) where a building is set
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'urban', 'suburban', 'exurban', 'rural', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**filter type** : ``filter_type``

* Definition: a device which removes solid particulates or airborne molecular contaminants
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'particulate air filter', 'chemical air filter', 'low-MERV pleated media', 'HEPA filter', 'electrostatic air treatment', 'gas-phase air treatment', 'ultraviolet air treatment', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**heating and cooling system type** : ``heat_cool_type``

* Definition: methods of conditioning or heating a room or building
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'radiant system', 'heat pump', 'forced air system', 'steam forced heat', 'wood stove', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**indoor space** : ``indoor_space``

* Definition: a distinguishable space within a structure, the purpose for which discrete areas of a building is used
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'bedroom', 'office', 'bathroom', 'foyer', 'kitchen', 'locker room', 'hallway', 'elevator', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**indoor surface** : ``indoor_surf``

* Definition: type of indoor surface
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'counter top', 'window', 'wall', 'cabinet', 'ceiling', 'door', 'shelving', 'vent cover', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**light type** : ``light_type``

* Definition: application of light to achieve some practical or aesthetic effect. Lighting includes the use of both artificial light sources such as lamps and light fixtures, as well as natural illumination by capturing daylight. Can also include absence of light
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'natural light', 'electric light', 'no light', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**occupancy at sampling** : ``occup_samp``

* Definition: number of occupants present at time of sample within the given space
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {integer}

**occupant density at sampling** : ``occupant_dens_samp``

* Definition: average number of occupants at time of sampling per square footage
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float}

**space typical state** : ``space_typ_state``

* Definition: customary or normal state of the space
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'typical occupied', 'typically unoccupied', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**substructure type** : ``substructure_type``

* Definition: the substructure or under building is that largely hidden section of the building which is built off the foundations to the ground floor level
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'crawlspace', 'slab on grade', 'basement', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**surface material** : ``surf_material``

* Definition: surface materials at the point of sampling
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'concrete', 'wood', 'stone', 'tile', 'plastic', 'glass', 'vinyl', 'metal', 'carpet', 'stainless steel', 'pint', 'cinder blocks', 'hay bales', 'stucco', 'adobe', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**surface-air contaminant** : ``surf_air_cont``

* Definition: contaminant identified on surface
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'dust', 'organic matter', 'particulate matter', 'volatile organic compounds', 'biological contaminants', 'radon', 'nutrients', 'biocides', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**typical occupant density** : ``typ_occupant_dens``

* Definition: customary or normal density of occupants
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float}

**ventilation rate** : ``ventilation_rate``

* Definition: ventilation rate of the system in the sampled premises
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: cubic meter per minute, liters per second

**ventilation type** : ``ventilation_type``

* Definition: ventilation system used in the sampled premises
* Field Format: free text
* Expected value: ventilation type name
* Value syntax: {text}
* Example: forced ventilation, mechanical ventilation, natural ventilation

link
~~~~

**link to classification information** : ``link_class_info``

* Definition: link to digitized soil maps or other soil classification information
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**link to climate information** : ``link_climate_info``

* Definition: link to climate resource
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**links to additional analysis** : ``link_addit_analys``

* Definition: link to additional analysis results performed on the sample
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

local environment conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**absolute air humidity** : ``abs_air_humidity``

* Definition: actual mass of water vapor - mh20 - present in the air water vapor mixture
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit} [kg|lb]
* Preferred unit: kilogram

**air temperature** : ``air_temp``

* Definition: temperature of the air at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit} [deg C]
* Preferred unit: degree Celsius

**alkalinity** : ``alkalinity``

* Definition: alkalinity, the ability of a solution to neutralize acids to the equivalence point of carbonate or bicarbonate
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milliequivalent per liter

**annual and seasonal precipitation** : ``annual_season_precpt``

* Definition: mean annual and seasonal precipitation (mm)
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millimeter

**annual and seasonal temperature** : ``annual_season_temp``

* Definition: mean annual and seasonal temperature (oC)
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree Celsius

**atmospheric data** : ``atmospheric_data``

* Definition: measurement of atmospheric data; can include multiple data
* Field Format: free text
* Expected value: atmospheric data name;measurement value
* Value syntax: {text};{float} {unit}

**barometric pressure** : ``barometric_press``

* Definition: force per unit area exerted against a surface by the weight of air above that surface
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: millibar

**climate environment** : ``climate_environment``

* Definition: treatment involving an exposure to a particular climate; can include multiple climates
* Field Format: free text
* Expected value: climate name;treatment duration;interval;experimental duration
* Value syntax: {text};{period};{interval};{period}

**conductivity** : ``conduc``

* Definition: electrical conductivity of water
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milliSiemens per centimeter

**current land use** : ``cur_land_use``

* Definition: present state of sample site
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'cities', 'farmstead', 'industrial areas', 'roads/railroads', 'rock', 'sand', 'gravel', 'mudflats', 'salt flats', 'badlands', 'permanent snow or ice', 'saline seeps', 'mines/quarries', 'oil waste areas', 'small grains', 'row crops', 'vegetable crops', 'horticultural plants (e.g. tulips)', 'marshlands (grass,sedges,rushes)', 'tundra (mosses,lichens)', 'rangeland', 'pastureland (grasslands used for livestock grazing)', 'hayland', 'meadows (grasses,alfalfa,fescue,bromegrass,timothy)', 'shrub land (e.g. mesquite,sage-brush,creosote bush,shrub oak,eucalyptus)', 'successional shrub land (tree saplings,hazels,sumacs,chokecherry,shrub dogwoods,blackberries)', 'shrub crops (blueberries,nursery ornamentals,filberts)', 'vine crops (grapes)', 'conifers (e.g. pine,spruce,fir,cypress)', 'hardwoods (e.g. oak,hickory,elm,aspen)', 'intermixed hardwood and conifers', 'tropical (e.g. mangrove,palms)', 'rainforest (evergreen forest receiving >406 cm annual rainfall)', 'swamp (permanent or semi-permanent water body dominated by woody plants)', 'crop trees (nuts,fruit,christmas trees,nursery trees)', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**current vegetation** : ``cur_vegetation``

* Definition: vegetation classification from one or more standard classification systems, or agricultural crop
* Field Format: free text
* Expected value: current vegetation type
* Value syntax: {text}

**current vegetation method** : ``cur_vegetation_meth``

* Definition: reference or method used in vegetation classification
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**dew point** : ``dew_point``

* Definition: the temperature to which a given parcel of humid air must be cooled, at constant barometric pressure, for water vapor to condense into water.
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree Celsius


**downward PAR** : ``down_par``

* Definition: visible waveband radiance and irradiance measurements in the water column
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: microEinstein per square meter per second

**drainage classification** : ``drainage_class``

* Definition: drainage classification from a standard system such as the USDA system
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'very poorly', 'poorly', 'somewhat poorly', 'moderately well', 'well', 'excessively drained', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**FAO classification** : ``fao_class``

* Definition: soil classification from the FAO World Reference Database for Soil Resources. The list can be found at http://www.fao.org/3/i3794en/I3794en.pdf
* Field Format: text choice
* Expected value: enumeration
* Value syntax: {term}

**fluorescence** : ``fluor``

* Definition: raw (volts) or converted (mg Chla/m^3) fluorescence of the water
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: milligram chlorophyll a per cubic meter

**gaseous environment** : ``gaseous_environment``

* Definition: use of conditions with differing gaseous environments; should include the name of gaseous compound, amount administered, treatment duration, interval and total experimental duration; can include multiple gaseous environment regimens
* Field Format: free text
* Expected value: gaseous compound name;gaseous compound amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: micromole per liter

**horizon** : ``horizon``

* Definition: specific layer in the land area which measures parallel to the soil surface and possesses physical characteristics which differ from the layers above and beneath
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'O horizon', 'A horizon', 'E horizon', 'B horizon', 'C horizon', 'R layer', 'Permafrost', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**horizon method** : ``horizon_meth``

* Definition: reference or method used in determining the horizon
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**humidity** : ``humidity``

* Definition: amount of water vapour in the air, at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: gram per cubic meter

**industrial effluent percent** : ``indust_eff_percent``

* Definition: percentage of industrial effluents received by wastewater treatment plant
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage

**light intensity** : ``light_intensity``

* Definition: measurement of light intensity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: lux

**local classification** : ``local_class``

* Definition: soil classification based on local soil classification system
* Field Format: free text
* Expected value: local classification name
* Value syntax: {text}

**local classification method** : ``local_class_meth``

* Definition: reference or method used in determining the local soil classification
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**mean friction velocity** : ``mean_frict_vel``

* Definition: measurement of mean friction velocity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: meter per second

**mean peak friction velocity** : ``mean_peak_frict_vel``

* Definition: measurement of mean peak friction velocity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: meter per second

**particle classification** : ``particle_class``

* Definition: particles are classified, based on their size, into six general categories:clay, silt, sand, gravel, cobbles, and boulders; should include amount of particle preceded by the name of the particle type; can include multiple values
* Field Format: free text
* Expected value: particle name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: micrometer

**pH** : ``ph``

* Definition: pH measurement
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float}

**pH method** : ``ph_meth``

* Definition: reference or method used in determining pH
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**photon flux** : ``photon_flux``

* Definition: measurement of photon flux
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: micromole per square meter per second

**pollutants** : ``pollutants``

* Definition: pollutant types and, amount or concentrations measured at the time of sampling; can report multiple pollutants by entering numeric values preceded by name of pollutant
* Field Format: free text
* Expected value: pollutant name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: gram, mole per liter, milligram per liter

**porosity** : ``porosity``

* Definition: porosity of deposited sediment is volume of voids divided by the total volume of sample
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage

**presence of pets or farm animals** : ``pet_farm_animal``

* Definition: specification of presence of pets or farm animals in the environment of subject, if yes the animals should be specified; can include multiple animals present
* Field Format: free text
* Expected value: presence status;type of animal or pet
* Value syntax: {boolean};{text}

**pressure** : ``pressure``

* Definition: pressure to which the sample is subject, in atmospheres
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: atmosphere

**relative air humidity** : ``rel_air_humidity``

* Definition: partial vapor and air pressure, density of the vapor and air, or by the actual mass of the vapor and air
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit} [%]
* Preferred unit: percentage

**seasonal environment** : ``season_environment``

* Definition: treatment involving an exposure to a particular season (e.g. winter, summer, rabi, rainy etc.)
* Field Format: free text
* Expected value: seasonal environment name;treatment duration;interval;experimental duration
* Value syntax: {text};{period};{interval};{period}

**sediment type** : ``sediment_type``

* Definition: information about the sediment type based on major constituents
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'biogenous', 'cosmogenous', 'hydrogenous', 'lithogenous', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**sewage type** : ``sewage_type``

* Definition: Type of sewage based on origin: wastewater treatment plant (municipal or industrial), open sewer line, river, stream, stagnant pool, or other.
* Field Format: free text
* Expected value: sewage type name
* Value syntax: {text}

**soil type** : ``soil_type``

* Definition: soil series name or other lower-level classification
* Field Format: free text
* Expected value: soil type name
* Value syntax: {text}

**soil type method** : ``soil_type_meth``

* Definition: reference or method used in determining soil series name or other lower-level classification
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**solar irradiance** : ``solar_irradiance``

* Definition: the amount of solar energy that arrives at a specific area of a surface during a specific time interval
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: watts per square meter, ergs per square centimeter per second

**surface humidity** : ``surf_humidity``

* Definition: surfaces: water activity as a function of air and material moisture
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit} [%]
* Preferred unit: percentage

**surface moisture** : ``surf_moisture``

* Definition: water held on a surface
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: parts per million, gram per cubic meter, gram per square meter

**surface moisture pH** : ``surf_moisture_ph``

* Definition: pH measurement of surface
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {integer [0-14]}

**surface temperature** : ``surf_temp``

* Definition: temperature of the surface at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit} [deg C]
* Preferred unit: degree Celsius

**temperature** : ``temp``

* Definition: temperature of the sample at time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree Celsius

**texture** : ``texture``

* Definition: the relative proportion of different grain sizes of mineral particles in a soil, as described using a standard system; express as % sand (50 um to 2 mm), silt (2 um to 50 um), and clay (<2 um) with textural name (e.g., silty clay loam) optional.
* Field Format: free text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: % sand/silt/clay

**texture method** : ``texture_meth``

* Definition: reference or method used in determining soil texture
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**tidal stage** : ``tidal_stage``

* Definition: stage of tide
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'low', 'high', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**total depth of water column** : ``tot_depth_water_col``

* Definition: measurement of total depth of water column
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: meter

**turbidity** : ``turbidity``

* Definition: turbidity measurement
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: formazin turbidity unit, formazin nephelometric units

**wastewater type** : ``wastewater_type``

* Definition: the origin of wastewater such as human waste, rainfall, storm drains, etc.
* Field Format: free text
* Expected value: wastewater type name
* Value syntax: {text}

**water current** : ``water_current``

* Definition: measurement of magnitude and direction of flow within a fluid
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: cubic meter per second, knots

**wind direction** : ``wind_direction``

* Definition: wind direction is the direction from which a wind originates
* Field Format: free text
* Expected value: wind direction name
* Value syntax: {text}

**wind speed** : ``wind_speed``

* Definition: speed of wind measured at the time of sampling
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: meter per second, kilometer per hour

local environment conditions imposed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**air temperature regimen** : ``air_temp_regm``

* Definition: information about treatment involving an exposure to varying temperatures; should include the temperature, treatment duration, interval and total experimental duration; can include different temperature regimens
* Field Format: free text
* Expected value: temperature value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: degree Celsius

**antibiotic regimen** : ``antibiotic_regm``

* Definition: information about treatment involving antibiotic administration; should include the name of antibiotic, amount administered, treatment duration, interval and total experimental duration; can include multiple antibiotic regimens
* Field Format: free text
* Expected value: antibiotic name;antibiotic amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: milligram

**chemical administration** : ``chem_administration``

* Definition: list of chemical compounds administered to the host or site where sampling occurred, and when (e.g. antibiotics, N fertilizer, air filter); can include multiple compounds. For Chemical Entities of Biological Interest ontology (CHEBI) (v111), please see http://purl.bioontology.org/ontology/CHEBI
* Field Format: free text
* Expected value: CHEBI;timestamp
* Value syntax: {term}; {timestamp}

**chemical mutagen** : ``chem_mutagen``

* Definition: treatment involving use of mutagens; should include the name of mutagen, amount administered, treatment duration, interval and total experimental duration; can include multiple mutagen regimens
* Field Format: free text
* Expected value: mutagen name;mutagen amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: milligram per liter

**fertilizer regimen** : ``fertilizer_regm``

* Definition: information about treatment involving the use of fertilizers; should include the name fertilizer, amount administered, treatment duration, interval and total experimental duration; can include multiple fertilizer regimens
* Field Format: free text
* Expected value: fertilizer name;fertilizer amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**fungicide regimen** : ``fungicide_regm``

* Definition: information about treatment involving use of fungicides; should include the name of fungicide, amount administered, treatment duration, interval and total experimental duration; can include multiple fungicide regimens
* Field Format: free text
* Expected value: fungicide name;fungicide amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**gravity** : ``gravity``

* Definition: information about treatment involving use of gravity factor to study various types of responses in presence, absence or modified levels of gravity; can include multiple treatments
* Field Format: free text
* Expected value: gravity factor value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: meter per square second, g

**growth hormone regimen** : ``growth_hormone_regm``

* Definition: information about treatment involving use of growth hormones; should include the name of growth hormone, amount administered, treatment duration, interval and total experimental duration; can include multiple growth hormone regimens
* Field Format: free text
* Expected value: growth hormone name;growth hormone amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**growth media** : ``growth_med``

* Definition: information about growth media for growing the plants or tissue cultured samples
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['','soil', 'liquid', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

**herbicide regimen** : ``herbicide_regm``

* Definition: information about treatment involving use of herbicides; information about treatment involving use of growth hormones; should include the name of herbicide, amount administered, treatment duration, interval and total experimental duration; can include multiple regimens
* Field Format: free text
* Expected value: herbicide name;herbicide amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**humidity regimen** : ``humidity_regm``

* Definition: information about treatment involving an exposure to varying degree of humidity; information about treatment involving use of growth hormones; should include amount of humidity administered, treatment duration, interval and total experimental duration; can include multiple regimens
* Field Format: free text
* Expected value: humidity value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: gram per cubic meter

**mechanical damage** : ``mechanical_damage``

* Definition: information about any mechanical damage exerted on the plant; can include multiple damages and sites
* Field Format: free text
* Expected value: damage type;body site
* Value syntax: {text};{text}

**mineral nutrient regimen** : ``mineral_nutr_regm``

* Definition: information about treatment involving the use of mineral supplements; should include the name of mineral nutrient, amount administered, treatment duration, interval and total experimental duration; can include multiple mineral nutrient regimens
* Field Format: free text
* Expected value: mineral nutrient name;mineral nutrient amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**non-mineral nutrient regimen** : ``non_mineral_nutr_regm``

* Definition: information about treatment involving the exposure of plant to non-mineral nutrient such as oxygen, hydrogen or carbon; should include the name of non-mineral nutrient, amount administered, treatment duration, interval and total experimental duration; can include multiple non-mineral nutrient regimens
* Field Format: free text
* Expected value: non-mineral nutrient name;non-mineral nutrient amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**perturbation** : ``perturbation``

* Definition: type of perturbation, e.g. chemical administration, physical disturbance, etc., coupled with time that perturbation occurred; can include multiple perturbation types
* Field Format: free text
* Expected value: perturbation type name;time interval
* Value syntax: {text};{interval}

**pesticide regimen** : ``pesticide_regm``

* Definition: information about treatment involving use of insecticides; should include the name of pesticide, amount administered, treatment duration, interval and total experimental duration; can include multiple pesticide regimens
* Field Format: free text
* Expected value: pesticide name;pesticide amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, mole per liter, milligram per liter

**pH regimen** : ``ph_regm``

* Definition: information about treatment involving exposure of plants to varying levels of pH of the growth media; can include multiple regimen
* Field Format: free text
* Expected value: measurement value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}

**radiation regimen** : ``radiation_regm``

* Definition: information about treatment involving exposure of plant or a plant part to a particular radiation regimen; should include the radiation type, amount or intensity administered, treatment duration, interval and total experimental duration; can include multiple radiation regimens
* Field Format: free text
* Expected value: radiation type name;radiation amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: becquerel

**rainfall regimen** : ``rainfall_regm``

* Definition: information about treatment involving an exposure to a given amount of rainfall; can include multiple regimens
* Field Format: free text
* Expected value: measurement value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: millimeter

**salt regimen** : ``salt_regm``

* Definition: information about treatment involving use of salts as supplement to liquid and soil growth media; should include the name of salt, amount administered, treatment duration, interval and total experimental duration; can include multiple salt regimens
* Field Format: free text
* Expected value: salt name;salt amount;treatment duration;interval;experimental duration
* Value syntax: {text};{float} {unit};{period};{interval};{period}
* Preferred unit: gram, microgram, mole per liter, gram per liter

**standing water regimen** : ``standing_water_regm``

* Definition: treatment involving an exposure to standing water during a plant's life span, types can be flood water or standing water; can include multiple regimens
* Field Format: free text
* Expected value: standing water type;treatment duration;interval;experimental duration
* Value syntax: {text};{period};{interval};{period}

**tissue culture growth media** : ``tiss_cult_growth_med``

* Definition: description of plant tissue culture growth media used
* Field Format: free text
* Expected value: PMID,DOI,url or free text
* Value syntax: {PMID\|DOI\|URL}

**water temperature regimen** : ``water_temp_regm``

* Definition: information about treatment involving an exposure to water with varying degree of temperature; can include multiple regimens
* Field Format: free text
* Expected value: measurement value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: degree Celsius

**watering regimen** : ``watering_regm``

* Definition: information about treatment involving an exposure to watering frequencies; can include multiple regimens
* Field Format: free text
* Expected value: measurement value;treatment duration;interval;experimental duration
* Value syntax: {float} {unit};{period};{interval};{period}
* Preferred unit: milliliter, liter

local environment history
~~~~~~~~~~~~~~~~~~~~~~~~~

**agrochemical additions** : ``agrochem_addition``

* Definition: addition of fertilizers, pesticides, etc. - amount and time of applications
* Field Format: free text
* Expected value: agrochemical name;agrochemical amount;timestamp
* Value syntax: {text};{float} {unit};{timestamp}
* Preferred unit: gram, mole per liter, milligram per liter

**crop rotation** : ``crop_rotation``

* Definition: whether or not crop is rotated, and if yes, rotation schedule
* Field Format: free text
* Expected value: crop rotation status;schedule
* Value syntax: {boolean};Rn/{timestamp}/{period}

**extreme event** : ``extreme_event``

* Definition: unusual physical events that may have affected microbial populations
* Field Format: free text
* Expected value: date
* Value syntax: {timestamp}

**fire** : ``fire``

* Definition: historical and/or physical evidence of fire
* Field Format: free text
* Expected value: date
* Value syntax: {timestamp}

**flooding** : ``flooding``

* Definition: historical and/or physical evidence of flooding
* Field Format: free text
* Expected value: date
* Value syntax: {timestamp}

**previous land use** : ``previous_land_use``

* Definition: previous land use and dates
* Field Format: free text
* Expected value: land use name;date
* Value syntax: {text};{timestamp}

**previous land use method** : ``previous_land_use_meth``

* Definition: reference or method used in determining previous land use and dates
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**tillage** : ``tillage``

* Definition: note method(s) used for tilling
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'drill', 'cutting disc', 'ridge till', 'strip tillage', 'zonal tillage', 'chisel', 'tined', 'mouldboard', 'disc plough', 'not applicable', 'not collected', 'not provided', 'restricted access', 'missing']

unusual properties
~~~~~~~~~~~~~~~~~~

**aluminium saturation** : ``al_sat``

* Definition: aluminum saturation (esp. for tropical soils)
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage

**aluminium saturation method** : ``al_sat_meth``

* Definition: reference or method used in determining Al saturation
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**extreme salinity** : ``extreme_salinity``

* Definition: measured salinity
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: practical salinity unit

**heavy metals** : ``heavy_metals``

* Definition: heavy metals present and concentrations of any drug used by subject and the frequency of usage; can include multiple heavy metals and concentrations
* Field Format: free text
* Expected value: heavy metal name;measurement value
* Value syntax: {text};{float} {unit}
* Preferred unit: microgram per gram

**heavy metals method** : ``heavy_metals_meth``

* Definition: reference or method used in determining heavy metals
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

**salinity method** : ``salinity_meth``

* Definition: reference or method used in determining salinity
* Field Format: free text
* Expected value: PMID,DOI or url
* Value syntax: {PMID\|DOI\|URL}

other
~~~~~

**description** : ``description``

* Definition: Description of the sample.
* Field Format: free text

**miscellaneous parameter** : ``misc_param``

* Definition: any other measurement performed or parameter collected, that is not listed here
* Field Format: free text
* Expected value: parameter name;measurement value
* Value syntax: {text};{float} {unit}


Standardised missing value vocabulary
-------------------------------------

not applicable
  information is inappropriate to report, can indicate that the standard itself fails to model or represent the information appropriately

not collected
  information of an expected format was not given because it has not been collected

not provided
  information of an expected format was not given, a value may be given at the later stage

restricted access
  information exists but can not be released openly because of privacy concerns

missing
  If you don't want to distinguish 'not collected', 'not provided' and 'restricted access'

Reference: https://press3.mcs.anl.gov/gensc/uncategorized/reporting-missing-values/


CNSA sample checklist
---------------------

- Pathogen affecting public health
   Use for pathogen samples that are relevant to public health.
  - Clinical or host-associated pathogen
  - Environmental, food or other pathogen
  - Combined pathogen: Batch submissions that include both clinical and environmental pathogen.
- Microbial sample
   Use for bacteria or other unicellular microbes when it is not appropriate or advantageous to use MIxS, Pathogen or Virus packages.
- Model organism or animal sample
   Use for multicellular samples or cell lines derived from common laboratory model organisms, e.g., mouse, rat, Drosophila, worm, fish, frog, or large mammals including zoo and farm animals.
- Metagenome or environmental sample
   Use for metagenomic and environmental samples when it is not appropriate or advantageous to use MIxS packages.
- Invertebrate sample
   Use for any invertebrate sample.
- Human sample
   Only use for human samples or cell lines that have no privacy concerns. For all studies involving human subjects, it is the submitter's responsibility to ensure that the information supplied protects participant privacy in accordance with all applicable laws, regulations and institutional policies. Make sure to remove any direct personal identifiers from your submission.

   For samples isolated from humans use the Pathogen, Microbe or appropriate MIxS package.
- Plant sample
   Use for any plant sample or cell line.
- Virus sample
   Use for all virus samples not directly associated with disease. Viral pathogens should be submitted using the Pathogen: Clinical or host-associated pathogen package.
- GSC MIxS environmental sample
   Genomic Standards Consortium package extension for reporting of measurements and observations obtained from the environment where the sample was obtained. The samples are validated for compliance based on the presence of the required core attributes as described in MIxS (http://gensc.org/gc_wiki/index.php/MIxS).
  - GSC MIxS air
  - GSC MIxS built environment
  - GSC MIxS host associated
  - GSC MIxS human associated
  - GSC MIxS human gut
  - GSC MIxS human oral
  - GSC MIxS human skin
  - GSC MIxS human vaginal
  - GCS MIxS microbial mat biolfilm
  - GSC MIxS miscellaneous natural or artificial environment
  - GSC MIxS plant associated
  - GSC MIxS sediment
  - GSC MIxS soil
  - GSC MIxS wastewater sludge
  - GSC MIxS water
- Other

  - Beta-lactamase
     Use for beta-lactamase gene transformants that have antibiotic resistance data.


Sample validation
-----------------

Typical validations include:

- Content must be supplied for mandatory fields. If information is unavailable for any mandatory field, please enter 'not applicable', 'not collected',  'not provided', 'restricted access' or 'missing' as appropriate, except for specific formats, please see field description for details.
- Multiple Samples cannot have identical attributes. You should have one Sample for each specimen, and each of your Samples must have differentiating information (excluding sample name, sample title, and description). This check was implemented to encourage submitters to include distinguishing information in their samples. If it is necessary to represent true biological replicates as separate Samples, you might add an 'aliquot' or 'replicate' attribute, e.g., 'replicate = biological replicate 1', as appropriate. Note that multiple assay types, e.g., RNA-seq and ChIP-seq data may reference the same Sample if appropriate.
- The values provided for some attributes are validated, for example, 'collection date' and 'geographic location' values must be provided in a recognized format. Refer to the attribute definition for information about required formats.

Organism validation rule

- Pathogen affecting public health: must have lineage Bacteria_, Viruses_ or Fungi_, to species-level
- Microbial sample: must have lineage Bacteria_, Archaea_, Viruses_, Viroids or Fungi_; or for unicellular eukaryotes, must not have lineage Metazoa_, Embryophyta_, `unclassified sequences`_ or `other sequences`_
- Model organism or animal sample: must NOT be Homo sapiens or have lineage Bacteria_, Archaea_, Viruses_ or Fungi_
- Metagenome or environmental sample: must be a metagenome, where lineage starts with unclassified sequences and scientific name ends with 'metagenome'
- Invertebrate sample: no check
- Human sample: must be Homo sapiens
- Plant sample: must have lineage Viridiplantae_, or be known to have plastids
- Virus sample: must have lineage Viruses_
- GSC MIxS environmental sample: no check
- Beta-lactamase: must have lineage Bacteria_

.. _Bacteria: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=2
.. _Archaea: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=2157
.. _Viruses: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=10239
.. _Fungi: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=4751
.. _Metazoa: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=33208
.. _Embryophyta: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=3193
.. _unclassified sequences: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=12908
.. _other sequences: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=28384
.. _Viridiplantae: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=33090
