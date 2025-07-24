# Clinical Trial Analysis: Brain Stimulation for Depression in Elderly Patients

This repository contains the full analysis of a randomized, controlled clinical trial evaluating the efficacy of non-invasive brain stimulation for the treatment of depression in elderly patients. 
The code is written in R.

## 📂 Project Structure

- (1) Preparation and cleaning of the dataset, including formatting, recoding, and generating long and wide data formats for later analyses.
- (2) Primary statistical analyses related to treatment effects on depression severity.
- (3) Analysis of secondary outcomes, including neuropsychological and psychological measures, as well as biomarkers.
- (0) Collection of custom functions used throughout the project, including utility functions, testing for overdispersion in GLMs and automated selection of Poisson or negative-binomial regression for count data.  

## 📊 Primary Analyses
- [**Rendered report (HTML)**](https://matthiasluthi.github.io/RCT-brain-stimulation-elderly/2-primary-analyses.html) - Full output with tables, figures, and statistical results

The primary statistical analyses include:
- Change in depression severity (e.g., HDRS-17) across time
- Response and remission rates
- Group comparisons using linear-mixed models with and without autoregressive covariance structure, regressions, GLMs, chi-square tests, etc.
- Analysis of blinding integrity
- Evaluation of adverse events associated with brain stimulation
- Various figures and tables, showing clinical and demographic data at baseline and development of clinical meausres during the study
- Automated creation and formatting of results tables and saving to Word documents, automated saving of figures

## 🧠 Secondary Analyses

Additional analyses were conducted on the following secondary outcomes:
- (a) Neuropsychology: Addenbrooke's Cognitive Examination–Revised (ACE-R) - [**Rendered report (HTML)**](https://matthiasluthi.github.io/RCT-brain-stimulation-elderly/3a-secondary-analyses-neuropsychology.html)
- (b) Psychology: Autobiographical Memory Test (AMT), Revised NEO Personality Inventory (NEO PI-R)- [**Rendered report (HTML)**](https://matthiasluthi.github.io/RCT-brain-stimulation-elderly/3b-secondary-analyses-psychology.html)
- (c) Biomarkers: Measures of motor cortical excitability related to the intervention

Each secondary outcome was tested for:
- Treatment-related change
- Associations with change in depression severity
- Predictive value of baseline measures for antidepressant response

## ⚠️ Disclaimer

Due to privacy and ethical restrictions, this repository does **not** include patient data.
Please contact the author for inquiries about collaboration or adaptation.
