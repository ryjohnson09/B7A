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

# Process Raw 16S reads using DADA2
# Depends on:	data/raw/microbiome_sequencing/ETEC/ETEC_sequences/F_reads/*.fastq.gz
#		data/raw/microbiome_sequencing/ETEC/ETEC_sequences/R_reads/*.fastq.gz
#		code/microbiome_analysis/DADA2_analysis_B7A.R
# Produces:	results/microbiome_analysis/F_reads_quality.png
#		results/microbiome_analysis/R_reads_quality.png
#		results/microbiome_analysis/filter_results.csv
#		results/microbiome_analysis/seqtab.rds
DADA2_objects =	results/microbiome_analysis/F_reads_quality.png\
		results/microbiome_analysis/R_reads_quality.png\
		results/microbiome_analysis/filter_results.csv\
		results/microbiome_analysis/seqtab.rds

all_DADA2: $(DADA2_objects)
.PHONY: all_DADA2

$(DADA2_objects) : data/raw/microbiome_sequencing/ETEC/ETEC_sequences/F_reads/*.fastq.gz\
		   data/raw/microbiome_sequencing/ETEC/ETEC_sequences/R_reads/*.fastq.gz\
		   code/microbiome_analysis/DADA2_analysis_B7A.R
	R -e "source('code/microbiome_analysis/DADA2_analysis_B7A.R', echo=T)"



# Length filter, chimera check, and assign taxonomy to DADA2 results
# Depends on:	data/processed/microbiome_analysis/silva_nr_v132_train_set.fa.gz
#		results/microbiome_analysis/seqtab.rds
#		code/microbiome_analysis/filter_chimera_taxonomy.R
# Produces:	results/microbiome_analysis/chimera_remaining.txt
#		results/microbiome_analysis/seqtab_final.rds
#		results/microbiome_analysis/tax_final.rds
filter_objects = results/microbiome_analysis/chimera_remaining.txt\
		 results/microbiome_analysis/seqtab_final.rds\
		 results/microbiome_analysis/tax_final.rds

all_filter: $(filter_objects)
.PHONY: all_filter

$(filter_objects) : data/processed/microbiome_analysis/silva_nr_v132_train_set.fa.gz\
                   results/microbiome_analysis/seqtab.rds\
                   code/microbiome_analysis/filter_chimera_taxonomy.R
	R -e "source('code/microbiome_analysis/filter_chimera_taxonomy.R', echo=T)"
