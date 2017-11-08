all.out <- list()
all.out[["Comment[ArrayExpressAccession]"]] <- "E-MTAB-5308"
all.out[["MAGE-TAB Version"]] <- "1.1"
all.out[["Investigation Title"]] <- "RNA-seq of human cell line HeLa after depletion of a lncRNA with three different loss-of-function methods"
all.out[["Comment[Submitted Name]"]] <- "RNA-seq of human cell line HeLa after depletion of a lncRNA with three different loss-of-function methods"
all.out[["Experiment Description"]] <- "Long noncoding RNAs (lncRNAs) are a major transcriptional output of the mammalian genome, and their cellular roles are typically assayed by a variety of loss-of-function approaches. This study aims to identify the best current method to deplete nuclear lncRNAs. Small interfering RNAs (RNAi), antisense oligonucleotides (LNAs) and CRISPR interference (CRISPRi) were applied to knock down loc100289019 (a typical nuclear lncRNA, referred to as lnc289) in HeLa cells. We generated sequencing libraries after performing each step of each method, up to and including depletion of lnc289. Differential expression analyses between libraries generated before and after each step allowed us to evaluate the effect of that step on gene expression. The transcriptional effect of lncRNA depletion was then compared to the magnitude of off-target effects inherent to each method."


all.out[["Experimental Design"]] <- c("genotype design",
                                      "compound treatment design",
                                      "quality control testing design",
                                      "cellular modification design")
all.out[["Experimental Design Term Source REF"]] <- c("EFO", 
                                                      "EFO",
                                                      "EFO",
                                                      "EFO")
all.out[["Experimental Design Term Accession Number"]] <- c("EFO_0001748",
                                                            "EFO_0001755",
                                                            "EFO_0001774",
                                                            "EFO_0004666")

all.out[["Experimental Factor Name"]] <- c("compound", 
                                           "genotype",
                                           "loss of function method")
all.out[["Experimental Factor Type"]] <- all.out[["Experimental Factor Name"]]
all.out[["Experimental Factor Term Source REF"]] <- ""
all.out[["Experimental Factor Term Accession Number"]] <- ""

all.out[["Person Last Name"]] <- "Lun"
all.out[["Person First Name"]] <- "Aaron"
all.out[["Person Mid Initials"]] <- "TL"
all.out[["Person Email"]] <- "aaron.lun@cruk.cam.ac.uk"
all.out[["Person Phone"]] <- ""
all.out[["Person Fax"]] <- ""      
all.out[["Person Address"]] <- "University of Cambridge Li Ka Shing Centre Robinson Way Cambridge CB2 0RE United Kingdom"
all.out[["Person Affiliation"]] <- "Cancer Research UK Cambridge Institute"
all.out[["Person Roles"]] <- "submitter"

all.out[["Protocol Name"]] <- c("P-MTAB-53234", 
                                "P-MTAB-53235",
                                "P-MTAB-53236",
                                "P-MTAB-53237",
                                "P-MTAB-53238",
                                "P-MTAB-53239",
                                "P-MTAB-53240",
                                "P-MTAB-53241")
all.out[["Protocol Type"]] <- c("sample collection protocol",
                                "growth protocol",
                                "treatment protocol",
                                "nucleic acid extraction protocol",
                                "nucleic acid library construction protocol",
                                "nucleic acid sequencing protocol",
                                "high throughput sequence alignment protocol",
                                "conversion protocol")
all.out[["Protocol Term Source REF"]] <- c("EFO", 
                                           "EFO",   
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO")
all.out[["Protocol Term Accession Number"]] <- c("EFO_0005518",
                                                 "EFO_0003789",
                                                 "EFO_0003969",
                                                 "EFO_0002944",
                                                 "EFO_0004184",
                                                 "EFO_0004170",
                                                 "EFO_0004917",
                                                 "EFO_0005520")

all.out[["Protocol Description"]] <- c("HeLa cells were obtained from American Type Culture Collection. To generate a heterogeneous population of cells expressing dCas9-KRAB, HeLa cells were transduced with a lentiviral vector containing the pHR-SFFV-dCAS9-BFP-KRAB vector along with polybrene (5 ug/ml, Sigma), incubated for 72 hours with medium replacement after 24 hours, and sorted for the BFP-expressing cells using a BD FACSAria III cell sorter. To generate single-cell clones expressing dCas9-KRAB, BFP-sorted single cells were plated on 96-well plates and grown for two weeks until colonies were obtained. Expression of dCas9-KRAB in each colony was verified by Western blot.", 
                                       "HeLa cells were maintained in Dulbecco's modified Eagle's medium (Sigma Aldrich, D6429) supplemented with 10% fetal bovine serum (Thermo Fisher Scientific) and cultured at 37 degrees Celsius with with 5% CO2.",
                                       "For RNAi and LNA, HeLa cells were transfected with Lipofectamine RNAiMax reagent (Thermo Fischer Scientific) following the manufacturer's instructions, using siRNAs (Thermo Fischer Scientifc) and LNA Gapmers (Exiqon) at a final concentration of 50 nM and 25nm, respectively. All experiments were done 48 hr after transfection. For CRISPRi, dCAS9-KRAB transduced cells were plated on a 12-well plate and infected with lentivirus containing the appropriate sgRNA vector. The virus was diluted with HeLa media (1:1 dilution) and cationic polymer polybrene was added to facilitate viral transduction (5 ug/ml, Sigma). After 24 hr incubation, supernatant was removed and fresh media was added for 48 hr before RNA collection. Samples with the same batch number were treated at roughly the same time (within the same month), while the experiment number is nested within batch and refers to cells treated on the same day.",
                                       "RNA (1 ug) was extracted with the RNeasy Kit (QIAGEN) and treated with DNase I following the manufacturer's instructions. RNA quality was assessed using a Total RNA Nano chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "RNA-seq libraries were prepared from HeLA cells using TruSeq Stranded Total RNA Kit with Ribo-Zero Gold (Illumina, RS-122-2303). Library quality was assessed using a DNA1000 chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "Indexed libraries were PCR amplified and sequenced on multiple lanes of an Illumina Hiseq 2500 instrument to obtain 125 bp paired-end reads.",
                                       "Reads were aligned to the hg38 build of the human genome using subread v1.5.3 in paired-end RNA-seq mode. The number of read pairs mapped to the exonic regions of each gene was then counted for each library, using the featureCounts function in Rsubread v1.28.0 with Ensembl GRCh38 version 90. Only alignments with mapping quality scores above 10 and forming reversely stranded fragments were considered during counting.",
                                       "The QuantiTect Reverse Transcription Kit (QIAGEN) was used for cDNA synthesis including an additional step to eliminate genomic DNA contamination.")
all.out[["Protocol Hardware"]] <- c("BD FACSAria III cell sorter",
                                    "", 
                                    "",
                                    "2100 Bioanalyzer",
                                    "2100 Bioanalyzer",
                                    "Illumina Hiseq 2500",
                                    "",
                                    "")
all.out[["Protocol Software"]] <- c("", 
                                    "",
                                    "",
                                    "",
                                    "",
                                    "",
                                    "(R)subread",
                                    "")
                                 
all.out[["Term Source Name"]] <- "EFO"
all.out[["Term Source File"]] <- "http://www.ebi.ac.uk/efo/"
all.out[["Term Source Version"]] <- ""
all.out[["Public Release Date"]] <- "2017-05-31"
all.out[["Comment[AEExperimentType]"]] <- "RNA-seq of coding RNA"
all.out[["Comment[SequenceDataURI]"]] <- "http://www.ebi.ac.uk/ena/data/view/ERR1751045-ERR1751697"
all.out[["Comment[SecondaryAccession]"]] <- "ERP020478"
all.out[["SDRF File"]] <- "E-MTAB-5308.sdrf.txt"

unlink("idf.tsv")
for (x in names(all.out)) {
    write(file="idf.tsv", paste0(c(x, all.out[[x]]), collapse="\t"), append=TRUE)
}


