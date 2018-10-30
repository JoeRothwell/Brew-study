# Brew study: analysis of commonly consumed coffee brews by metabolomics
## Exploratory analysis of untargeted data
## Targeted analysis of known coffee compounds
Data extracted from raw data files with ProFinder. A custom PCDL was used as a target for the extraction of 120 coffee compounds found in online databases. Once the data was extracted, each compound group was checked for poorly drawn peaks. Bad peaks were redrawn or removed if judged to be noise. Blanks were checked to make sure extracted peaks were not spurious. Any judged to be noise only were deleted.
Selected samples, rejected samples and blanks were extracted (182 samples)
Pooled QCs were extracted separately for calculation of CVs.
88 Compound groups were extracted from the coffee compounds initially collected. 17 were rejected
#for the following reasons, leaving 71.
*1 due to detector saturation: caffeine
*2 due to not of interest: Ochratoxin, 2 isomers
*2 due to high single ion compounds: 2-Ac-5-Methylthiophene, 2-Me-5-(1-propenyl)pyrazine
*8 due to high levels in blanks: Kahweol, 3-(3,4-dihydroxyphenyl)-2-propenoic acid, Linoleic acid, 
*4-Ethyl-1,2-dimethoxybenzene, Campestanol, Trimethylamine, Feruloylquinic acid (3)
*1 due to too many missing peaks: Liberine
*1 due to mistaken identity: Trigonelline (RT is around 0.7 min)
*3 due to high CV in pooled QCs: diacetin, Nicotinic acid, Phenylalanine
### Update at June 2018: More compounds excluded
File was exported from ProFinder as "export detailed csv". Due to the commas in compound names, it wasndifficult to parse the csv to be read in directly. Therefore, the csv was manually manipulated in Excel to correct the shifted rows. Codes were given to each extracted compound eg cmpd_152_3.11. Excel file location: \\Inti\BMA\Coffee project NCI\Targeted analysis 88 cmpds Feb 15.xlsx.The "cleaned" sheet was exported to the R project directory.
