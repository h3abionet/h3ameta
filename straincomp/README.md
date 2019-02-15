# h3ameta - straincomp

## Challenge:
You have sequence data generated from a culture of bacteria isolated from the bloodstream of an individual with an infection. You also have sequence data from their stool.
1. How similar are the bacteria in their bloodstream and their stool?
2. Identify nucleotide differences between sequences found in the isolate data and the bloodstream data, and determine whether the sequences found in the isolate could have originated$
3. Identify antibiotic resistance genes in the genome of the bloodstream isolate and in the gut metagenomes. (Bonus)

## 1- Development Plan

### 1. Processes
#### a- [Classification](https://github.com/h3abionet/h3ameta/tree/master/examples/taxonomic_classification) (Taxonomic)
#### b- [Strainsifter](https://github.com/tamburinif/StrainSifter) (Relatedness)
#### c- [SRST2](https://github.com/katholt/srst2) (Antimicrobial Resistance)

![classification_srst2_flowchart"](https://github.com/h3abionet/h3ameta/blob/master/straincomp/classification_srst2_dag.png "classification_srst2_flowchart")

### 2. Implementation assigments
~~- to implement **SRST2** in the **classification.nt**~~ (Done :+1:)
- to implement **Strainsifter** in the **classification_srst2.nt** :computer: :sweat:

### 3. Timeline
~~- by **14 Feb 2019** should have the first draft of the pipeline **classificaon_srst2.nt**~~ (Done :+1:)
- by **01 Mar 2019** should have the complete draft of the pipeline **classificaon_srst2_strainsifter.nt** :computer: :sweat:

### 4. Communication Plan
Slack

### 5. Members
- **Heyam Mohammed** (University of Khartoum, Sudan) 👩🏻 🇸🇩
- **Mushal Allam** (National Institute for Communicable Diseases, South Africa) 👨🏽 🇿🇦
- **Penistacia Maela** (University of Johannesburg, South Africa) 👩🏻 🇿🇦
