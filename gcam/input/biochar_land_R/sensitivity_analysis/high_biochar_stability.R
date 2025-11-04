devtools::load_all()

# supply sector   subsector       technology
# biochar         slow pyrolysis  poultry biochar
# biochar         slow pyrolysis  dairy biochar
# biochar         slow pyrolysis  beef biochar
# biochar         slow pyrolysis  pork biochar
# biochar         slow pyrolysis  goat biochar

# conversions
# 1e6 tons = 1 Mt
# 1e12 MJ = 1 EJ

### ADDITIONAL FILES DATA + SOURCES
## A_an_secout
#               source 1: Bentley, J. A. et al. Economics of Dairy Manure Management in Iowa. Iowa State University Animal Industry Report 13, (2016).
#               source 2: Malone, G. W. Nutrient Enrichment in integrated Broiler Production Systems. Poultry Science 71, 1117–1122 (1992).
#               source 3: https://www.fao.org/3/i6421e/i6421e.pdf
#               source 5: https://www.fao.org/3/t0279e/T0279E05.htm - guidelines for meat cutting and processing - table 3
#               source 9: Lefebvre, D. et al. Biomass residue to carbon dioxide removal: quantifying the global impact of manure_fuel. manure_fuel 5, 65 (2023).
#        poultry (units): 0.0273 (kg collectible manure solids /day /head) [9] * 1 (head) / 1.008 (kg meat/head) [3] * 51 day lifespan [2] = 1.381 Mt manure/Mt Poultry
#          dairy (units): 0.5391 (kg collectible manure solids /day /head) [9] * 1 (head) * 365 (days) / (23,578 (lbs milk/head/year) [1] * 0.453592 (kg/lbs))= 0.018 Mt manure/Mt dairy product
#           beef (units): 0.5391 (kg collectible manure solids /day /head) [9] * 1 (head) / 228 (kg meat/head) [5] * 3 (year lifespan) * 365 (days) [5] = 2.589 Mt manure/Mt Beef
#           pork (units): 0.0930 (kg collectible manure solids /day /head) [9] * 1 (head) / 55.76 (kg meat/head) [5] * 0.5 (year lifespan) * 365 (days) [5] = 0.304 Mt manure/Mt Pork
#           goat (units): 0.2433 (kg collectible manure solids /day /head) [9] * 1 (head) / 19.1 (kg meat/head) [5] * .667 (year lifespan) * 365 (days) [5] = 3.101 Mt manure/Mt Goat
# calculations in plant_costs.xlsx

## A_An_secout_prices:
# manure prices are set to 1e-4, as that is the smallest possible price based on the rounding in L202.an_input

#A21.globaltech_cost.csv is assumed to model capital costs, per description in A22.globaltech_cost_low.csv, given that the same columns are present
# instead of adding additional costs on the biochar input (pre-treatment/post-treatment costs for biochar are considered as part of the capital costs) - per docs /supply_energy
example_file <- find_csv_file("energy/A21.globaltech_cost", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)

#           source 1: plant costs.xlsx
tmp[15] <- "biochar,slow pyrolysis,poultry_biochar,non-energy,0.0251,0.0251"
tmp[16] <- "biochar,slow pyrolysis,pork_biochar,non-energy,0.0252,0.0252"
tmp[17] <- "biochar,slow pyrolysis,beef_biochar,non-energy,0.0254,0.0254"
tmp[18] <- "biochar,slow pyrolysis,dairy_biochar,non-energy,0.0274,0.0274"
tmp[19] <- "biochar,slow pyrolysis,goat_biochar,non-energy,0.0242,0.0242"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



#modify A21.globaltech_coef.csv
example_file <- find_csv_file("energy/A21.globaltech_coef", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)

#               source 1: Baniasadi, M. et al. Waste to energy valorization of poultry litter by slow pyrolysis. Renewable Energy 90, 458–468 (2016).
#               source 2: Zhou, S., Liang, H., Han, L., Huang, G. and Yang, Z., 2019. The influence of manure feedstock, slow pyrolysis, and hydrothermal temperature on manure thermochemical and combustion properties. Waste management, 88, pp.85-95.
#               source 4: Zeng, X., Xiao, Z., Zhang, G., Wang, A., Li, Z., Liu, Y., Wang, H., Zeng, Q., Liang, Y. and Zou, D., 2018. Speciation and bioavailability of heavy metals in pyrolytic biochar of swine and goat manures. Journal of Analytical and Applied Pyrolysis, 132, pp.82-93.
# original value (units): unitless
#     used value (units): 47% yield implies 2.1276 Mt feedstock / Mt biochar
tmp[16] <- "biochar,slow pyrolysis,poultry_biochar,poultry manure,2.139,2.139,2.139" #[1]
tmp[17] <- "biochar,slow pyrolysis,pork_biochar,pork manure,2.136,2.136,2.136" # [4]
tmp[18] <- "biochar,slow pyrolysis,beef_biochar,beef manure,2.150,2.150,2.150" #[2]
tmp[19] <- "biochar,slow pyrolysis,dairy_biochar,dairy manure,2.325,2.325,2.325" #[2]
tmp[20] <- "biochar,slow pyrolysis,goat_biochar,goat manure,2.055,2.055,2.055" #[4]
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)


#modify A21.globaltech_shrwt.csv
example_file <- find_csv_file("energy/A21.globaltech_shrwt", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
#               source 1: A322.globaltech_shrwt
# original value (units): same assumptions as fertilizer production from hydrogen
#     used value (units): same assumptions as fertilizer production from hydrogen
tmp[15] <- "biochar,slow pyrolysis,poultry_biochar,1,1"
tmp[16] <- "biochar,slow pyrolysis,pork_biochar,1,1"
tmp[17] <- "biochar,slow pyrolysis,beef_biochar,1,1"
tmp[18] <- "biochar,slow pyrolysis,dairy_biochar,1,1"
tmp[19] <- "biochar,slow pyrolysis,goat_biochar,1,1"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



#modify A21.globaltech_interp.csv
example_file <- find_csv_file("energy/A21.globaltech_interp", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)

#               source 1: A322.globaltech_shrwt
#               source 2: Woolf, D., Amonette, J. E., Street-Perrott, F. A., Lehmann, J. & Joseph, S. Sustainable biochar to mitigate global climate change. Nat Commun 1, 56 (2010).
# original value (units): s-curve following [2]
#     used value (units): same assumptions as fertilizer production from hydrogen
tmp[9]  <- "biochar,slow pyrolysis,poultry_biochar,share-weight,final-calibration-year,end-year,,fixed"
tmp[10] <- "biochar,slow pyrolysis,pork_biochar,share-weight,final-calibration-year,end-year,,fixed"
tmp[11] <- "biochar,slow pyrolysis,beef_biochar,share-weight,final-calibration-year,end-year,,fixed"
tmp[12] <- "biochar,slow pyrolysis,dairy_biochar,share-weight,final-calibration-year,end-year,,fixed"
tmp[13] <- "biochar,slow pyrolysis,goat_biochar,share-weight,final-calibration-year,end-year,,fixed"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



#modify A21.subsector_interp
example_file <- find_csv_file("energy/A21.subsector_interp", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
#               source 1: A322.globaltech_shrwt
# original value (units): same assumptions as fertilizer production from hydrogen
#     used value (units): same assumptions as fertilizer production from hydrogen
tmp[11] <- "biochar,slow pyrolysis,share-weight,final-calibration-year,end-year,,fixed,0"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



#modify A21.subsector_shrwt
example_file <- find_csv_file("energy/A21.subsector_shrwt", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
#               source 1: A322.globaltech_shrwt
# original value (units): same assumptions as fertilizer production from hydrogen
#     used value (units): same assumptions as fertilizer production from hydrogen
tmp[12] <- "biochar,slow pyrolysis,final-calibration-year,,1,0"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



#modify A21.sector
example_file <- find_csv_file("energy/A21.sector", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)

#               source 1: A21.Sector
# original value (units): Mt, because the primary output of this sector is the biochar which is used as a mass input to fertilizer production, and input from poultry is also in Mt
#     used value (units): otherwise, use same logit.exponents as the other considered technologies.
tmp[13] <- "biochar,Mt,Mt,1975$/kg,-3,0," #primary output is biochar
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)


#modify A21.subsector_logit
example_file <- find_csv_file("energy/A21.subsector_logit", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
#               source 1: A21.subsector_logit
# original value (units): same values as other supply sectors
#     used value (units):
tmp[12] <- "biochar,slow pyrolysis,-6,0,absolute-cost-logit"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)



### CHANGES TO EMISSIONS ###
#modify A_PrimaryFuelCCoef.csv
example_file <- find_csv_file("emissions/A_PrimaryFuelCCoef.csv", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
#     used value (units): C stored as biochar
#               source 1: Woolf, D., Amonette, J. E., Street-Perrott, F. A., Lehmann, J. & Joseph, S. Sustainable biochar to mitigate global climate change. Nat Commun 1, 56 (2010).

### THESE VALUES ARE FOR CARBON STORED IN BIOCHAR ###
### values calculated in 41467_2010 excel spreadsheet
tmp[37] <- "beef manure,-0.0393,0"
tmp[38] <- "dairy manure,-0.0393,0"
tmp[39] <- "goat manure,-0.448,0"
tmp[40] <- "pork manure,-0.448,0"
tmp[41] <- "poultry manure,-0.458,0"
print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)


#modify A41.tech_coeff.csv
example_file <- find_csv_file("emissions/A41.tech_coeff.csv", FALSE)[[1]]
file.copy(from = example_file, to = paste0(example_file, ".bak"))

# Change one value in file, then rewrite to same path
tmp <- readr::read_lines(example_file)
print("\n file before changes")
print(tmp)
# original value (units): (Mt CH4 per Mt biochar supply)
#               source 1: Woolf, D., Amonette, J. E., Street-Perrott, F. A., Lehmann, J. & Joseph, S. Sustainable biochar to mitigate global climate change. Nat Commun 1, 56 (2010).

### THESE VALUES ARE FOR CARBON AVOIDANCE AND ARE NOT SUBJECT TO CARBON SUBSIDIES ###
### values calculated in 41467_2010 excel spreadsheet
tmp[41] <- "biochar,slow pyrolysis,dairy_biochar,,,,-0.0105,-0.00052,,,,,,"
tmp[42] <- "biochar,slow pyrolysis,beef_biochar,,,,-0.0105,-0.00052,,,,,,"
tmp[43] <- "biochar,slow pyrolysis,pork_biochar,,,,-0.1023,-0.00061,,,,,,"
tmp[44] <- "biochar,slow pyrolysis,goat_biochar,,,,-0.1023,-0.00061,,,,,,"
tmp[45] <- "biochar,slow pyrolysis,poultry_biochar,,,,-0.0057,-0.00074,,,,,,"

print("\n file after changes")
print(tmp)
readr::write_lines(tmp, example_file)

### update xml files ###
devtools::load_all()
driver_drake()
