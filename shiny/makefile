###############################################
#### List all targes by typing `make list` ####
###############################################

.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs -n 1





########################
#### Tidy Data Sets ####
########################


# Create Clinical Metadata Table Extracted from TrEAT DB
# Depends on:	data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX
#               data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX
#               code/Create_Clin_Metadata.R
# Produces:     data/processed/TrEAT_Clinical_Metadata_tidy.csv
data/processed/TrEAT_Clinical_Metadata_tidy.csv : data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX\
                                                  data/raw/TrEAT_Merge_DataDictionary_2018.06.27.XLSX\
                                                  code/Create_Clin_Metadata.R
	R -e "source('code/Create_Clin_Metadata.R', echo=T)"


# Convert raw Olink data to tidy format
# Depends on:	data/raw/20170276_Henry_M_Jackson_Foundation-Ventura_NPX_LOD_Updated_and_Revised_2.26.18.xlsx
#		code/Process_Raw_Olink.R
# Produces:	data/processed/Olink_tidy.csv
data/processed/Olink_tidy.csv : data/raw/20170276_Henry_M_Jackson_Foundation-Ventura_NPX_LOD_Updated_and_Revised_2.26.18.xlsx\
				code/Process_Raw_Olink.R
	R -e "source('code/Process_Raw_Olink.R', echo=T)"
