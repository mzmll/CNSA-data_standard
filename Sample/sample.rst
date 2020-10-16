Sample
========

Description of biological source material; each physically unique specimen should be registered as a single sample with a unique set of attributes.

Sample classification
---------------------

geography
~~~~~~~~~

**profile position(profile_position)**

* Definition: cross-sectional position in the hillslope where sample was collected. sample area position in relation to surrounding areas
* Field Format: text choice
* Expected value: enumeration
* Value syntax: ['', 'summit', 'shoulder', 'backslope', 'footslope', 'toeslope']


**slope aspect(slope_aspect)**

* Definition: the direction a slope faces. While looking down a slope use a compass to record the direction you are facing (direction or degrees); e.g., NW or 315 degrees. This measure provides an indication of sun and wind exposure that will influence soil temperature and evapotranspiration.
* Field Format: free text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: degree

**slope gradient(slope_gradient)**

* Definition: commonly called 'slope'. The angle between ground surface and a horizontal line (in percent). This is the direction that overland water would flow. This measure is usually taken with a hand level meter or clinometer
* Field Format: restricted text
* Expected value: measurement value
* Value syntax: {float} {unit}
* Preferred unit: percentage


non-sample terms
~~~~~~~~~~~~~~~~

**project accession(project_accession)**

* Definition: A valid project accession has 'CNP' prefix.
* Field Format: restricted text
* Expected value: CNSA project accession

**experimental factor(experimental_factor)**

* Definition: Experimental factors are essentially the variable aspects of an experiment design which can be used to describe an experiment, or set of experiments, in an increasingly detailed manner. This field accepts ontology terms from Experimental Factor Ontology (EFO) and/or Ontology for Biomedical Investigations (OBI). For a browser of EFO (v 2.95) terms, please see http://purl.bioontology.org/ontology/EFO; for a browser of OBI (v 2018-02-12) terms please see http://purl.bioontology.org/ontology/OBI
* Field Format: free text
* Expected value: text or EFO and/or OBI
* Value syntax: {termLabel} {[termID]}|{text}
* Example: time series design [EFO:EFO_0001779]

**ploidy(ploidy)**

* Definition: The ploidy level of the genome (e.g. allopolyploid, haploid, diploid, triploid, tetraploid). It has implications for the downstream study of duplicated gene and regions of the genomes (and perhaps for difficulties in assembly). For terms, please select terms listed under class ploidy (PATO:001374) of Phenotypic Quality Ontology (PATO), and for a browser of PATO (v1.269) please refer to http://purl.bioontology.org/ontology/PATO
* Field Format: free text
* Expected value: PATO
* Value syntax: {term}
* Example: allopolyploid, polyploid

**number of replicons(num_replicons)**

* Definition: Reports the number of replicons in a nuclear genome of eukaryotes, in the genome of a bacterium or archaea or the number of segments in a segmented virus. Always applied to the haploid chromosome count of a eukaryote.
* Field Format: restricted text
* Expected value: for eukaryotes and bacteria: chromosomes (haploid count); for viruses: segments
* Value syntax: {integer}

**extrachromosomal elements(extrachrom_elements)**

* Definition: Do plasmids exist of significant phenotypic consequence (e.g. ones that determine virulence or antibiotic resistance). Megaplasmids? Other plasmids (borrelia has 15+ plasmids)
* Field Format: restricted text
* Expected value: number of extrachromosmal elements
* Value syntax: {integer}

**estimated size(estimated_size)**

* Definition: The estimated size of the genome (in bp) prior to sequencing. Of particular importance in the sequencing of (eukaryotic) genome which could remain in draft form for a long or unspecified period.
* Field Format: restricted text
* Expected value: number of base pairs
* Value syntax: {integer} bp
* Example: 300000 bp

****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:



****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:



****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:



****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:



****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:



****

* Definition:
* Field Format:
* Expected value:
* Value syntax:
* Preferred unit:
