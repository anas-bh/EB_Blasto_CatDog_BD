# Statistical Analysis — Molecular Epidemiology of E. bieneusi and Blastocystis in Cats and Dogs

## Associated publication
Harun AB, Bayazid AA, Hasan MF, et al. Molecular Epidemiology of Enterocytozoon bieneusi
and Blastocystis in Cats and Dogs in Dhaka and Gazipur, Bangladesh.
Food and Waterborne Parasitology (under review).

## Repository contents
- `EB_Blasto_Cat_Dog_AnasBH.R` — R script for all statistical analyses
- `CatDog_AnasBH.xlsx` — Epidemiological dataset (Sheet 1: Cat, Sheet 2: Dog)

## How to run
1. Install R and RStudio
2. Place both files in the same folder
3. Open the R script and run from top to bottom
4. Required packages: readxl, logistf, pROC, randomForest

## Analysis includes
- Chi-square univariable analysis
- Firth penalized likelihood logistic regression (Heinze & Schemper, 2002)
- ROC curve analysis (Supplementary Figures S1–S2)
- Random Forest variable importance (Supplementary Figures S1–S2)
