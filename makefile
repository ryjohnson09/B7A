###############################################
#### List all targes by typing `make list` ####
###############################################

.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs -n 1





################################
### Processing Raw Data Sets ###
################################

# Merge and Clean Raw Humichip Data
# Depends on:	data/raw/Humichip/HumiChipBx123-LTO.txt
#		data/raw/Humichip/HumiChipBxRpt-LTO.txt
#		code/Merge_Humichip_B7A.R
# Produces:	data/processed/Merged_humichip_B7A.tsv
data/processed/Merged_humichip_B7A.tsv : data/raw/Humichip/HumiChipBx123-LTO.txt\
			                 data/raw/Humichip/HumiChipBxRpt-LTO.txt\
			                 code/Merge_Humichip_B7A.R
	R -e "source('code/Merge_Humichip_B7A.R', echo=T)"


# Merge and Clean Raw Geochip Data
# Depends on:   data/raw/Geochip/GeoChipBx1-LTO.txt
#               data/raw/Geochip/GeoChipBx2-LTO.txt
#		data/raw/Geochip/GeoChipBx3-LTO.txt
#		data/raw/Geochip/GeoChipBxRpt-LTO.txt
#               code/Merge_Geochip_B7A.R
# Produces:     data/processed/Merged_geochip_B7A.tsv
data/processed/Merged_geochip_B7A.tsv : data/raw/Geochip/GeoChipBx1-LTO.txt\
			                data/raw/Geochip/GeoChipBx2-LTO.txt\
               				data/raw/Geochip/GeoChipBx3-LTO.txt\
			                data/raw/Geochip/GeoChipBxRpt-LTO.txt\
			                code/Merge_Geochip_B7A.R
	R -e "source('code/Merge_Geochip_B7A.R', echo=T)"




# Format B7A metadata
# Depends on:	data/raw/B7A_Listing_(Glomics)_10OCT2018.xlsx
#               data/processed/Merged_humichip_B7A.tsv
#		data/processed/Merged_geochip_B7A.tsv
#               code/Process_B7A_Metadata.R
# Produces:     data/processed/B7A_metadata.tsv
data/processed/B7A_metadata.tsv : data/raw/B7A_Listing_(Glomics)_10OCT2018.xlsx\
                                  data/processed/Merged_humichip_B7A.tsv\
				  data/processed/Merged_geochip_B7A.tsv\
                                  code/Process_B7A_Metadata.R
	R -e "source('code/Process_B7A_Metadata.R', echo=T)"








###########################
### Microbiome Analysis ###
###########################

# Raw data directory
raw_dir = data/raw/microbiome_sequencing/
# Processed data directory
processed_dir = data/processed/microbiome_analysis/
# Code
code_dir = code/microbiome_analysis/
# Results
results_dir = results/microbiome_analysis/


# Process Raw 16S reads using DADA2
# Depends on:	$(raw_dir)ETEC/ETEC_sequences/F_reads/*.fastq.gz
#		$(raw_dir)ETEC/ETEC_sequences/R_reads/*.fastq.gz
#		$(code_dir)DADA2_analysis_B7A.R
# Produces:	$(results_dir)F_reads_quality.png
#		$(results_dir)R_reads_quality.png
#		$(results_dir)filter_results.csv
#		$(results_dir)seqtab.rds
DADA2_objects =	$(results_dir)F_reads_quality.png\
		$(results_dir)R_reads_quality.png\
		$(results_dir)filter_results.csv\
		$(results_dir)seqtab.rds

all_DADA2: $(DADA2_objects)
.PHONY: all_DADA2

$(DADA2_objects) : $(raw_dir)ETEC/ETEC_sequences/F_reads/*.fastq.gz\
		   $(raw_dir)ETEC/ETEC_sequences/R_reads/*.fastq.gz\
		   $(code_dir)DADA2_analysis_B7A.R
	R -e "source('$(code_dir)DADA2_analysis_B7A.R', echo=T)"



# Length filter, chimera check, and assign taxonomy to DADA2 results
# Depends on:	$(processed_dir)silva_nr_v132_train_set.fa.gz
#		$(results_dir)seqtab.rds
#		$(code_dir)filter_chimera_taxonomy.R
# Produces:	$(results_dir)chimera_remaining.txt
#		$(results_dir)seqtab_final.rds
#		$(results_dir)tax_final.rds
filter_objects = $(results_dir)chimera_remaining.txt\
		 $(results_dir)seqtab_final.rds\
		 $(results_dir)tax_final.rds

all_filter: $(filter_objects)
.PHONY: all_filter

$(filter_objects) : $(processed_dir)silva_nr_v132_train_set.fa.gz\
                    $(results_dir)seqtab.rds\
                    $(code_dir)filter_chimera_taxonomy.R
	R -e "source('$(code_dir)filter_chimera_taxonomy.R', echo=T)"


# Generate Metadata to be used for phyloseq
# Depends on:	$(raw_dir)ED_B7A01_Microbiome_MergedWithClinicalData_2019.01.11.xlsx
#		$(raw_dir)MissingSamples_ClinicalData_2019.10.07.xlsx
#		$(raw_dir)B7A AbxTimes_09OCT2019.xlsx
#		$(code_dir)metadata_creation.R
# Produces:	$(processed_dir)B7A_metadata_clean.csv
$(processed_dir)B7A_metadata_clean.csv : $(raw_dir)ED_B7A01_Microbiome_MergedWithClinicalData_2019.01.11.xlsx\
        				 $(raw_dir)MissingSamples_ClinicalData_2019.10.07.xlsx\
					 $(raw_dir)B7A\ AbxTimes_09OCT2019.xlsx\
				    	 $(code_dir)metadata_creation.R
	R -e "source('$(code_dir)metadata_creation.R', echo=T)"


# Filter DADA2 data using Phyloseq
# Depends on:	$(results_dir)seqtab_final.rds
#               $(results_dir)tax_final.rds
#		$(processed_dir)B7A_metadata_clean.csv
#		$(code_dir)phyloseq_filtering.R
# Produces: 	$(results_dir)prev_plot_before_filt.png
#		$(results_dir)ps2.rds
$(results_dir)ps2.rds $(results_dir)prev_plot_before_filt.png : $(results_dir)seqtab_final.rds\
								$(results_dir)tax_final.rds\
								$(processed_dir)B7A_metadata_clean.csv\
								$(code_dir)phyloseq_filtering.R
	R -e "source('$(code_dir)phyloseq_filtering.R', echo=T)"


# Abundance Plots
# Depends on:	$(results_dir)ps2.rds
# 		$(code_dir)abundance_barplots.R
# Produces:	$(barplot_dir)barplot_phylum_patient.png
#               $(barplot_dir)barplot_phylum_time.png
#               $(barplot_dir)barplot_phylum_mod_sev.png
#               $(barplot_dir)barplot_phylum_naive.png
#               $(barplot_dir)barplot_phylum_abx.png
#               $(barplot_dir)barplot_genus_patient.png
#               $(barplot_dir)barplot_genus_mod_sev.png
#               $(barplot_dir)barplot_genus_time.png
#               $(barplot_dir)barplot_genus_naive.png
barplot_dir = results/microbiome_analysis/barplots/

barplots = $(barplot_dir)barplot_phylum_patient.png\
           $(barplot_dir)barplot_phylum_time.png\
           $(barplot_dir)barplot_phylum_mod_sev.png\
           $(barplot_dir)barplot_phylum_naive.png\
           $(barplot_dir)barplot_phylum_abx.png\
           $(barplot_dir)barplot_genus_patient.png\
           $(barplot_dir)barplot_genus_mod_sev.png\
           $(barplot_dir)barplot_genus_time.png\
           $(barplot_dir)barplot_genus_naive.png

all_barplots: $(barplots)
.PHONY: all_barplots

$(barplots) : $(results_dir)ps2.rds\
	      $(code_dir)abundance_barplots.R
	R -e "source('$(code_dir)abundance_barplots.R', echo=T)"
