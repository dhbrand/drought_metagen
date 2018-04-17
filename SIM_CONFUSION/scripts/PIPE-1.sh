#!/bin/bash

mkdir -p QC0
mkdir -p QC1
mkdir -p FILTER0
mkdir -p FILTER1
mkdir -p JOINs
mkdir -p COMBINE
mkdir -p FILTER2
mkdir -p HUMANN2
mkdir -p PATHOSCOPE2

fastqs="SampleAmmonia_assimilation SampleArginine_Biosynthesis_extended SampleBacterial_Chemotaxis SampleBacterial_Cytoskeleton SampleBiotin_biosynthesis SamplecAMP_signaling_in_bacteria SampleCobalamin_synthesis SampleCobalt-zinc-cadmium_resistance SampleCoenzyme_A_Biosynthesis SampleCoenzyme_B12_biosynthesis SampleCopper_homeostasis SampleCysteine_Biosynthesis SampleDe_Novo_Purine_Biosynthesis SampleDe_Novo_Pyrimidine_Synthesis SampleDNA_repair._bacterial SampleDNA_Repair_Base_Excision SampleDNA_repair._UvrABC_system SampleDNA-replication SampleECF_class_transporters SampleEntner-Doudoroff_Pathway"

for s in ${fastqs}
do
    sbatch -o ${s}.pipe.%A.out scripts/pipe.sh ${s}
done


