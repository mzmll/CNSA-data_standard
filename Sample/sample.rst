Sample
========

Description of biological source material; each physically unique specimen should be registered as a single sample with a unique set of attributes.

Sample classification
---------------------

geography
~~~~~~~~~

**profile position** : ``profile_position``

* Definition: cross-sectional position in the hillslope where sample was collected. sample area position in relation to surrounding areas
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'summit', 'shoulder', 'backslope', 'footslope', 'toeslope']


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

internal environment
~~~~~~~~~~~~~~~~~~~~

**building occupancy type** : ``build_occup_type``

* Definition: the primary function for which a building or discrete part of a building is intended to be used
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'airport', 'agricultural', 'commercial', 'educational', 'government', 'health care', 'high rise', 'industrial', 'low rise', 'market', 'office', 'military', 'parking', 'residential', 'restaurant', 'school', 'sports complex', 'storage', 'religious', 'transport', 'wood framed', 'missing', 'not applicable', 'not collected']

**building setting** : ``building_setting``

* Definition: a location (geography) where a building is set
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'urban', 'suburban', 'exurban', 'rural', 'missing', 'not applicable', 'not collected']

**filter type** : ``filter_type``

* Definition: a device which removes solid particulates or airborne molecular contaminants
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'particulate air filter', 'chemical air filter', 'low-MERV pleated media', 'HEPA filter', 'electrostatic air treatment', 'gas-phase air treatment', 'ultraviolet air treatment', 'missing', 'not applicable', 'not collected']

**heating and cooling system type** : ``heat_cool_type``

* Definition: methods of conditioning or heating a room or building
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'radiant system', 'heat pump', 'forced air system', 'steam forced heat', 'wood stove', 'missing', 'not applicable', 'not collected']

**indoor space** : ``indoor_space``

* Definition: a distinguishable space within a structure, the purpose for which discrete areas of a building is used
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'bedroom', 'office', 'bathroom', 'foyer', 'kitchen', 'locker room', 'hallway', 'elevator', 'missing', 'not applicable', 'not collected']

**indoor surface** : ``indoor_surf``

* Definition: type of indoor surface
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'counter top', 'window', 'wall', 'cabinet', 'ceiling', 'door', 'shelving', 'vent cover']

**light type** : ``light_type``

* Definition: application of light to achieve some practical or aesthetic effect. Lighting includes the use of both artificial light sources such as lamps and light fixtures, as well as natural illumination by capturing daylight. Can also include absence of light
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'natural light', 'electric light', 'no light', 'missing', 'not applicable', 'not collected']

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
* Value syntax: ['', 'typical occupied', 'typically unoccupied', 'missing', 'not applicable', 'not collected']

**substructure type** : ``substructure_type``

* Definition: the substructure or under building is that largely hidden section of the building which is built off the foundations to the ground floor level
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'crawlspace', 'slab on grade', 'basement']

**surface material** : ``surf_material``

* Definition: surface materials at the point of sampling
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'concrete', 'wood', 'stone', 'tile', 'plastic', 'glass', 'vinyl', 'metal', 'carpet', 'stainless steel', 'paint', 'cinder blocks', 'hay bales', 'stucco', 'adobe']

**surface-air contaminant** : ``surf_air_cont``

* Definition: contaminant identified on surface
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'dust', 'organic matter', 'particulate matter', 'volatile organic compounds', 'biological contaminants', 'radon', 'nutrients', 'biocides']

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
* Value syntax: ['', 'aerobic', 'anaerobic']

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
* Value syntax: ['', 'metagenome (not viral targeted)', 'viral fraction metagenome (virome)', 'sequence-targeted metagenome', 'metatranscriptome (not viral targeted)', 'viral fraction RNA metagenome (RNA virome)', 'sequence-targeted RNA metagenome', 'microbial single amplified genome (SAG)', 'viral single amplified genome (vSAG)', 'isolate microbial genome', 'missing', 'not applicable', 'not collected']
* Example: viral fraction metagenome (virome)

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
* Value syntax: ['', 'clean catch', 'catheter']

**virus enrichment approach** : ``virus_enrich_appr``

* Definition: Approach used to enrich the sample for viruses, if any. If more than one approach was used, include multiple ‘virus_enrich_appr’ fields.
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'filtration', 'ultrafiltration', 'centrifugation', 'ultracentrifugation', 'PEG Precipitation', 'FeCl Precipitation', 'CsCl density gradient', 'DNAse', 'RNAse', 'targeted sequence capture', 'missing', 'not applicable', 'not collected']
* Example: filtration + FeCl Precipitation + ultracentrifugation + DNAse

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

host description
~~~~~~~~~~~~~~~~

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:


**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:


**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:

**** : ````

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:
