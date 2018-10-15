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


