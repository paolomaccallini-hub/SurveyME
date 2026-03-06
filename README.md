# Reanalysis of the EMEA Pan-European ME Patient Survey by Machine Learning, Bimodal analysis, and Local Outlier Factor

## Abstract

Raw data from the EMEA Pan-European ME Patient Survey were analysed by several methods, including local outlier factor (LOF), supervised and unsupervised Machine Learning (ML), and bimodal fit of age at first symptoms. A total of 9,596 subjects were considered, including three patients added to the original data set. Age at first symptoms for males with ME/CFS displays a statistically significant fit with a gamma bimodal density, with about 60% of the subjects developing the disease at a mean age of 27 (sd = 13) and 40% at a mean age of 40 (sd = 8). Males who develop the disease earlier in their lives experience more sensitivity, dizziness, and sleep problems, but less pain. They show a higher level of severity. Females report more comorbidities than males, they have more PEM, and are less likely to be stable. Females are also less likely to begin the disease after an infection, but more likely to indicate vaccination as the triggering event. Despite these differences, females are not more severe and more fatigued than males. LOF analysis identified 480 outliers (5%) who display an increased prevalence of deterioration over time, are less likely to develop the disease after an infection, and have first symptoms at a younger age. Logistic regression and Naive Bayes models were built for the identification of housebound patients (AUC 84%), outliers (AUC 88%), and subjects who tend to deteriorate over time (AUC 76%). The application of both LOF and Naive Bayes models trained on it shows promise for the identification of atypical ME patients who are candidate subjects for the investigation of monogenic causes of this disabling condition. 

## Methods

### Data source

The results of the 2021 survey promoted by the European ME Alliance (EMEA) were kindly shared with me by one of the authors, upon request, along with a copy of the form. Since the raw data were never made public by EMEA, I won't share them in this repository. I received a `.xlsx` file including the answers of 11,297 patients. A description of the survey and a detailed analysis of the results have been released by EMEA in April 2024 ([Angelsen A et Schei T, 2024](https://www.europeanmealliance.org/emea-pan-european-survey-UK.shtml)).

### Data filtering and editing

The original raw data included 11,297 subjects and 272 variables. Individuals without a ME/CFS diagnosis were removed (1,484 subjects), then 69 subjects with no specified gender were filtered out. Then, 4 patients whose year of birth was after the year of first symptoms were removed. Two males indicated as trigger pregnancy/birth, but their illness started years after birth: they were removed. Of the remaining 9,738 subjects, two male patients reported endometrosis as a comorbidity, and their data were edited to remove the inconsistency. Missing data were present in 145 subjects, and they were removed. The subjects selected from the original data set are 9,593.

Age at diagnosis was calculated as the year of first symptoms minus the year of birth. Age was calculated as 2021 minus the year of birth. Three further cases collected by me were added to the raw data (folder [Supplementary Cases](Supplementary%20Cases/)), for a total of 9,596 ME/CFS patients. After editing and selecting the variables, I included 23 parameters for each subject. Below is the description of the variables of the present analysis.

| Variable                          | Type      | Values |
|:----------------------------------|:----------|:--------------------------------------------------------------------------------------------------------------------------------|
| gender_n                          | binary    | male (1), female (2) |
| trigger_n                         | integer   | infection (1), accident/injury/surgery (2), traumatic life event (3), vaccine (4), pregnancy/birth (5), none (6), other (7) |
| relatives_n                       | binary    | first degree relatives with ME/CFS (1), no relatives with ME/CFS (0) |
| fatigue_n                         | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| muscleorjointpain_n               | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| dizzness_n                        | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| sleepproblems_n                   | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| sensitivity_n                     | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| pem_n                             | ordinal   | not at all (1), a little (2), moderately (3), a lot (4), very much (5), not relevant (6) |
| course_n                          | ordinal   | mostly stable, small fluctuations (1); major fluctuations (2); fluctuating initially, then mostly stable (3); fluctuating initially, then mostly improvement (4); fluctuating initially, then mostly deterioration (5); mainly improvement (6); mainly deterioration (7) |
| fibromyalgia_n                    | binary    | yes (1), no (0) |
| Hashimotosthyroiditis_n           | binary    | yes (1), no (0) |
| interstitialcystitis_n            | binary    | yes (1), no (0) |
| siccasyndrome_n                   | binary    | yes (1), no (0) |
| irritablebowelsyndrome_n          | binary    | yes (1), no (0) |
| migraine_n                        | binary    | yes (1), no (0) |
| allergies_n                       | binary    | yes (1), no (0) |
| endometriosis_n                   | binary    | yes (1), no (0) |
| asthma_n                          | binary    | yes (1), no (0) |
| multiplechemicalsensitivities_n   | binary    | yes (1), no (0) |
| severity_n                        | ordinal   | completely recovered (1); better than mild (2); mild (~50% reduction) (3); moderate, mostly housebound (4); severe, mostly bedridden (5); very severe, totally bedridden, needs help (6) |
| age_ill                           | numeric   | age at illness onset (year_ill − year_born) |
| age                               | numeric   | age in years (2021 − year_born) |
| duration                          | numeric   | disease duration in years (2021 − year_ill) |

<p align="left">
  <em>Table 1. Variables of the present study. </em>
</p>

For the majority of the following analyses, I expanded variables `course_n` and `trigger_n` into several binary variables, as indicated in the following table.

| Binary Variable  | Original variable      | Original values |
|:-----------------|:-----------------------|:----------------|
|improveming_n     | course_n               | 4, 6            |
|deteriorating_n   | course_n               | 5, 7            |
|stable_n          | course_n               | 1, 3            |
|fluctuating_n     | course_n               | 2               |
|improveming_n     | course_n               | 4, 6            |
|infection_n       | trigger_n              | 1               |
|AIS_n             | trigger_n              | 2               |
|trauma_n          | trigger_n              | 3               |
|vaccine_n         | trigger_n              | 4               |
|pregnancy_n       | trigger_n              | 5               |
|no_trigger_n      | trigger_n              | 6               |
|other_n           | trigger_n              | 7               |

### Bimodal analysis of age at onset

The distributions of `age_ill` for females, males, and both sexes were analysed by maximum likelihood estimation (MLE) performed by the function `mixfit()` of the R package `mixR` ([You Y, 2021](https://joss.theoj.org/papers/10.21105/joss.04031)). I asked for a bimodal density using both gamma and lognormal densities. Then, I tested the fit by the two-sided Kolmogorov-Smirnov test ([Marsaglia G et al. 2003](https://www.jstatsoft.org/article/view/v008i18)), using the function `ks.test()` of the package `stats`. This method tests the following null hypothesis: the distribution comes from the fitted density. Therefore, a high p-value supports the fit.

### Local outlier factor

I searched for outliers among patients by using the local outlier factor (LOF) algorithm ([Breunig MM, 2000](https://dl.acm.org/doi/10.1145/335191.335388)). I set 10 as the minimum value for minPts (it represents the minimum size of a cluster for the determination of local outliers), and I chose 5% of the total number of patients as the maximum value, considering that it represents the size of the larger cluster that can be regarded as a cluster of outliers. I assigned an LOF for each patient, for each value of minPts included in this interval. For each patient, the maximum LOF was assigned as definitve value. The calculation of LOF was performed by the function `lof()` of the package `dbscan`. The distribution of the LOF was used to test the hypothesis of being an outlier for each of the three added cases by one-sided testing with a cut-off for significance of 0.05.

### Correlations

Correlation tables between variables were computed using the Spearman coefficient. Only significant correlations were reported (two-sided test, cut-off for significance of 0.05). I employed the package `corplot`. P

### Comparisons

Subjects were divided into two groups by various criteria (males vs females, outliers vs non-outliers, etc), then the two groups were compared using the Fisher test for binary variables and the T-test for non-binary ones. P-values were corrected for multiple comparisons using the Benjamini-Hochberg method, implemented by the function `p.adjust` of the package `stats`. 

### Supervised classification

I trained both a Naive Bayes and a Logistic regression model to recognise individuals who are mostly housebound to very severe (`severity_n` above 4). For both models, I estimated ean sensitivity, specificity, and accuracy as a function of the threshold of the predicted probability used for classification. I used the same approach to test logistic regression and Naive Bayes to predict patients who deteriorate. I then applied the two logistic models to the three added patients.

### Unsupervised classification and PCA

To search for clusters in an unbiased way, I used the function `pam()` of the package `cluster`. Also, Principal Component Analysis (PCA) was performed by the function `princomp()` of the package `stats`, and visualisation of the principal planes and distribution of the variance across the components by both a custom script and the  function `biplot()` of the same package.

## Results

### Bimodal analysis of age at onset

The bimodal fit based on two gamma densities for the age at onset was significant only for males (Figure 1), according to the two-sided Kolmogorov-Smirnov test. The lognormal fit gave a similar result, but the fit is less significant (see repository files).  

<img src="https://github.com/user-attachments/assets/02c41852-73d0-4de7-8853-3f6e46ccfe58" width="700" alt="Age_at_onset_gamma_Bimodality">
<p align="left">
  <em>Figure 1. Gamma bimodal fit of age at first symptoms, performed by maximum likelihood estimation (MLE) by the function `mixfit()` of the R package `mixR`. The two-sided Kolmogorov-Smirnov test was used to test the null: the distribution comes from the fitted density. </em>
</p>

The parameters of the two gamma densities for females, males, nd both sexes are reported in Table 2. For the analogous results of the lognormal distribution, see the files in the repository. Remember that only for males, the bimodal density significantly fits the data.

| Group               | Proportion | Mean  | SD   |
|---------------------|-----------:|------:|-----:|
| Younger Females     | 0.53       | 26.44 | 12.6 |
| Older Females       | 0.47       | 41.28 | 8.26 |
| Younger Males       | 0.62       | 26.71 | 12.77 |
| Older Males         | 0.38       | 42.75 | 8.15 |
| Younger Both sexes  | 0.54       | 26.44 | 12.61 |
| Older Both sexes    | 0.46       | 41.43 | 8.29 |

<p align="left">
  <em>Table 2. Parameters of the two gamma densities whose combination fits the distribution of age at first symptoms. Stratification by sex. </em>
</p>

### Local outlier factor and PCA

The distribution of the Local Outlier Factor (LOF) is shown in Figure 2, where the LOF of patient PT1 (folder [Supplementary Cases](Supplementary%20Cases/)) is also indicated by a red line. We see that this patient is an outlier because it displays an LOF that is greater than the 0.95 quantile of the distribution. This patient is a 46 yo male who got sick in his twenties and who has been housebound for most of the disease course. He also has low myalgia and high dizziness. According to this analysis, either he has another disease, or he is a case of monogenic ME/CFS (see discussion section). For the analogous plots of the other two added patients, see the file section.

<img src="https://github.com/user-attachments/assets/315558bb-ed4f-4ac9-adf6-c2cc6cd79477" width="1000" alt="LOF_histogram_9554">
<p align="left">
  <em>Figure 2. Distribution of Local Outlier Factor for the about 9,600 patients included in this analysis. The vertical line indicates the LOF of PT1 (one of the three added patients), a male subject who developed ME/CFS at 20 and who has been mostly housebound ever since. The p-value refers to the null: PT1 belongs to the population of ME/CFS patients. </em>
</p>

For each added patient, the script plots the position on the three principal planes after computing loadings (available in the file section). In Figure 3, the results of the PCA, with the position of PT1 on the principal planes.

<img src="https://github.com/user-attachments/assets/616f9871-9fe2-41aa-a491-358b2adce1c6" width="1000" alt="PCA PT1">
<p align="left">
  <em>Figure 3. Results of PCA, represented twice. Left: positions of all the patients on the three principal planes, including PT1 in red, and also representation of the space PC1-PC2-PC3. Right: same as before, with the projections of the axes of the space of the data on the principal planes (in red), the projection of PT1 (in red), and the distribution of the variance across the principal components. Loadings are available among the files of this repository. Since the variance is spread across most of the principal components, this analysis has little significance. </em>
</p>

### Correlations and comparisons

The correlation table is reported in Figure 4. Note that the two multi-level variables `course_n` and `trigger_n` were translated into binary variables, one for each level. Only significant correlations are reported (significance level: 0.05). We note that fatigue, PEM, pain, dizziness, sleep problems, sensitivity, and fibromyalgia positively correlate with deterioration. We note that females tend to correlate with all the comorbidities and most symptoms, except fatigue. 

<img src="https://github.com/user-attachments/assets/5dbacd61-fefc-456a-a45a-7ff5ebd385e3" width="1000" alt="EMEA_course_triggers_corr">
<p align="left">
  <em>Figure 4. Correlation table. Gender is 2 for females and 1 for males. AIS_n: accident/injury/surgery. </em>
</p>

When comparing males and females, the variables that are significantly different after correction for multiple comparisons are in Table 3. Females have more comorbidities overall, they are older, with a longer disease duration, have more PEM, and are less likely to be stable. Females are also less likely to begin the disease after an infection, but more likely to indicate vaccination as the triggering event (perhaps because of papillomavirus vaccination?). Females have a smaller LOF than males. Despite these differences, females are not more severe and more fatigued than males. The prevalence of endometriosis among women with ME/CFS does not appear to be higher than what has been documented in the general population ([Moradi Y et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9131783/))

| Variable | Females (n=8198) | Males (n=1398) | P-value | Adjusted P-value |
|----------|-----------------|----------------|---------|-----------------|
| endometriosis_n | 8.94% | 0% | 8.37e-53 | 1.42e-51 |
| sensitivity_n | 3.66 (1.1) | 3.19 (1.24) | 5.02e-39 | 5.69e-38 |
| fibromyalgia_n | 33.9% | 19.6% | 5.14e-28 | 4.37e-27 |
| migraine_n | 31.5% | 18.53% | 2.77e-24 | 1.88e-23 |
| pregnancy_n | 2.82% | 0.07% | 7.29e-15 | 4.13e-14 |
| muscleorjointpain_n | 3.62 (1.16) | 3.34 (1.31) | 2.07e-13 | 1.01e-12 |
| Hashimotosthyroiditis_n | 9.56% | 4.22% | 2.37e-12 | 1.01e-11 |
| irritablebowelsyndrome_n | 43.1% | 33.26% | 3.91e-12 | 1.48e-11 |
| interstitialcystitis_n | 3.5% | 0.93% | 9.57e-09 | 3.25e-08 |
| dizzness_n | 3.02 (1.25) | 2.81 (1.26) | 2.87e-08 | 8.87e-08 |
| siccasyndrome_n | 5.57% | 2.58% | 5.29e-07 | 1.50e-06 |
| age | 47.9 (13.43) | 45.99 (14.98) | 8.10e-06 | 2.12e-05 |
| duration | 15.49 (10.99) | 14.15 (10.64) | 1.46e-05 | 3.55e-05 |
| allergies_n | 38.55% | 32.55% | 1.69e-05 | 3.83e-05 |
| lof | 1.08 (0.12) | 1.1 (0.14) | 5.99e-05 | 1.27e-04 |
| multiplechemicalsensitivities_n | 14.2% | 10.44% | 1.23e-04 | 2.46e-04 |
| infection_n | 49.23% | 54.58% | 2.36e-04 | 4.46e-04 |
| vaccine_n | 4.44% | 2.43% | 2.71e-04 | 4.85e-04 |
| sleepproblems_n | 3.86 (1.09) | 3.75 (1.15) | 6.97e-04 | 1.18e-03 |
| trauma_n | 5.61% | 3.51% | 7.72e-04 | 1.25e-03 |
| pem_n | 4.62 (0.72) | 4.54 (0.81) | 1.60e-03 | 2.47e-03 |
| asthma_n | 16.59% | 13.45% | 3.15e-03 | 4.66e-03 |
| no_trigger_n | 13.88% | 16.88% | 3.80e-03 | 5.17e-03 |
| other_n | 13.88% | 16.88% | 3.80e-03 | 5.17e-03 |
| AIS_n | 4.7% | 3.22% | 1.18e-02 | 1.54e-02 |
| stable_n | 23.32% | 25.89% | 3.79e-02 | 4.77e-02 |

<p align="left">
  <em> Table 3. Comparison between females and males. Only significant differences were included, after BH correction for multiple comparisons (see file section for the complete table). Note that pregnancy_n indicates as trigger both pregnancy and birth. </em>
</p>

When we assign males to two groups based on the probability of belonging to the two densities in Table 1, we obtain the results collected in Table 4. We see that those who got sick at a younger age have a longer disease duration despite being younger. They are more severe and more isolated in the space of the variables, according to their mean LOF. They have more sensitivity, dizziness, and sleep problems, but less pain.

| Variable              | M_comp_1 (n=769)   | M_comp_2 (n=629)   | P-value      | Adjusted P-value |
|:----------------------|-------------------:|-------------------:|-------------:|-----------------:|
| age_ill               |  22.13 (8.84)      |  43.71 (6.57)      | 0.00e+00     | 0.00e+00         |
| comp.prob.1           |   0.88 (0.16)      |   0.30 (0.07)      | 0.00e+00     | 0.00e+00         |
| comp.prob.2           |   0.12 (0.16)      |   0.70 (0.07)      | 0.00e+00     | 0.00e+00         |
| age                   |  38.05 (14.27)     |  55.69 (8.89)      | 1.79e-137    | 1.52e-136        |
| duration              |  15.92 (12.13)     |  11.98 (7.97)      | 5.72e-13     | 3.89e-12         |
| severity_n            |   3.99 (0.87)      |   3.77 (0.81)      | 8.66e-07     | 4.91e-06         |
| lof                   |   1.11 (0.16)      |   1.09 (0.10)      | 1.20e-03     | 5.83e-03         |
| sensitivity_n         |   3.27 (1.25)      |   3.08 (1.22)      | 4.89e-03     | 2.08e-02         |
| dizzness_n            |   2.89 (1.28)      |   2.72 (1.22)      | 8.54e-03     | 2.90e-02         |
| sleepproblems_n       |   3.82 (1.13)      |   3.66 (1.17)      | 7.84e-03     | 2.90e-02         |
| muscleorjointpain_n   |   3.26 (1.35)      |   3.44 (1.27)      | 1.09e-02     | 3.37e-02         |

<p align="left">
  <em> Table 4. Comparison between males more likely to belong to those who developed the disease at a young age, versus the others, according to the results in Figure 1 and Table 2. Variable comp.prob.1 indicates the probability of being part of the density with a lower mean. Analogously for comp.prob.2. </em>
</p>

In Table 5, I collected the results of the comparison between the subjects with a severity from 4 (mostly housebound) to 6 (very severe), and all the others. The first group displais worse fatigue, dizziness, PEM, sensitivity, pain, fibromyalgia, sleep problems, migraine, allergies, asthma, interstitial cystitis, endometriosis, and sicca syndrome. The more severe group is more likely to be deteriorating, is younger, includes more females, developed the disease at a younger age, and has a longer disease duration. Also, patients with a more severe ME/CFS are more likely to have relatives with ME/CFS. Infection is more likely to be the initial trigger in housebound patients than in the less severe ones.

| Variable                      | Housebound (n=7093) | notHousebound (n=2503) | P-value      | Adjusted P-value |
|:------------------------------|--------------------:|-----------------------:|-------------:|-----------------:|
| severity_n                    |      4.29 (0.53)    |           2.86 (0.37)  | 0.00e+00     | 0.00e+00         |
| fatigue_n                     |      4.71 (0.59)    |           3.94 (0.85)  | 8.48e-313    | 1.44e-311        |
| deteriorating_n               |     55.67 %         |          17.66 %       | 9.21e-254    | 1.04e-252        |
| dizzness_n                    |      3.18 (1.24)    |           2.42 (1.10)  | 4.73e-170    | 4.02e-169        |
| pem_n                         |      4.75 (0.61)    |           4.20 (0.90)  | 9.15e-155    | 6.22e-154        |
| sensitivity_n                 |      3.77 (1.07)    |           3.07 (1.15)  | 1.29e-144    | 7.31e-144        |
| improving_n                   |      3.10 %         |          18.58 %       | 2.75e-126    | 1.34e-125        |
| muscleorjointpain_n           |      3.73 (1.16)    |           3.15 (1.16)  | 2.72e-97     | 1.16e-96         |
| sleepproblems_n               |      3.98 (1.04)    |           3.46 (1.16)  | 1.68e-84     | 6.35e-84         |
| stable_n                      |     19.70 %         |          35.04 %       | 1.60e-51     | 5.44e-51         |
| fibromyalgia_n                |     34.34 %         |          24.65 %       | 1.11e-19     | 3.43e-19         |
| irritablebowelsyndrome_n      |     43.96 %         |          35.16 %       | 1.18e-14     | 3.34e-14         |
| migraine_n                    |     31.58 %         |          24.01 %       | 5.39e-13     | 1.41e-12         |
| fluctuating_n                 |     21.53 %         |          28.73 %       | 5.84e-13     | 1.42e-12         |
| multiplechemicalsensitivities_n |   15.10 %         |           9.55 %       | 7.82e-13     | 1.77e-12         |
| age_ill                       |     31.88 (12.90)  |          33.58 (12.52) | 8.18e-09     | 1.74e-08         |
| lof                           |      1.08 (0.12)    |           1.10 (0.12)  | 2.14e-07     | 4.28e-07         |
| allergies_n                   |     38.80 %         |          34.48 %       | 1.24e-04     | 2.34e-04         |
| duration                      |     15.54 (11.11)  |          14.61 (10.43) | 1.64e-04     | 2.93e-04         |
| relatives_n                   |     15.56 %         |          13.14 %       | 3.33e-03     | 5.66e-03         |
| asthma_n                      |     16.73 %         |          14.42 %       | 6.57e-03     | 1.06e-02         |
| no_trigger_n                  |     13.76 %         |          15.90 %       | 9.61e-03     | 1.42e-02         |
| other_n                       |     13.76 %         |          15.90 %       | 9.61e-03     | 1.42e-02         |
| trauma_n                      |      4.96 %         |           6.27 %       | 1.28e-02     | 1.81e-02         |
| gender_n                      |     85.96 %         |          83.94 %       | 1.47e-02     | 2.00e-02         |
| age                           |     47.42 (13.61)  |          48.19 (13.89) | 1.75e-02     | 2.29e-02         |
| infection_n                   |     50.71 %         |          48.02 %       | 2.14e-02     | 2.69e-02         |
| interstitialcystitis_n        |      3.37 %         |           2.44 %       | 2.29e-02     | 2.78e-02         |
| endometriosis_n               |      7.99 %         |           6.63 %       | 2.86e-02     | 3.35e-02         |
| siccasyndrome_n               |      5.43 %         |           4.31 %       | 3.08e-02     | 3.49e-02         |
| AIS_n                         |      4.74 %         |           3.76 %       | 4.30e-02     | 4.72e-02         |

<p align="left">
  <em> Table 5. Comparison between individuals who are housebound or worse, and the others. AIS_n: accident/injury/surgery. </em>
</p>

ME/CFS patients with affected relatives (Table 6) have a longer disease duration and developed the disease earlier. They have more MCS, fibromyalgia, allergies, migraine, sensitivity, and dizziness. They are more severe and are more isolated in the variable space (higher LOF). They are older.

| Variables                        | aff_relatives (1433) | no_aff_relatives (8163) | p.val     | p.adjust  |
|----------------------------------|---------------------:|------------------------:|----------:|----------:|
| relatives_n                      | 100 %                | 0 %                     | 0.00e+00  | 0.00e+00  |
| duration                         | 18.2 (12.36)         | 14.78 (10.6)            | 2.39e-22  | 4.06e-21  |
| age_ill                          | 30.29 (13.73)        | 32.68 (12.63)           | 8.42e-10  | 9.54e-09  |
| lof                              | 1.11 (0.14)          | 1.08 (0.12)             | 1.36e-09  | 1.16e-08  |
| multiplechemicalsensitivities_n  | 18.63 %              | 12.78 %                 | 9.79e-09  | 6.66e-08  |
| fibromyalgia_n                   | 36.64 %              | 30.97 %                 | 2.85e-05  | 1.62e-04  |
| allergies_n                      | 41.87 %              | 36.93 %                 | 4.34e-04  | 2.11e-03  |
| migraine_n                       | 33.43 %              | 28.94 %                 | 7.00e-04  | 2.97e-03  |
| sensitivity_n                    | 3.68 (1.12)          | 3.57 (1.13)             | 9.57e-04  | 3.62e-03  |
| severity_n                       | 3.98 (0.82)          | 3.91 (0.79)             | 2.74e-03  | 9.32e-03  |
| dizzness_n                       | 3.07 (1.23)          | 2.97 (1.25)             | 6.15e-03  | 1.90e-02  |
| age                              | 48.49 (14.48)        | 47.47 (13.54)           | 1.31e-02  | 3.71e-02  |

<p align="left">
  <em> Table 6. Comparison between patients with and patients without affected relatives. </em>
</p>

In Table 7, I compared outliers (LOF above the 0.9 quantile) with the others. Outliers display an increased prevalence of deterioration over time, are less likely to develop the disease after an infection, and have first symptoms at a younger age. They have less fatigue and less PEM. Variances are larger among outliers, suggesting that they do not define a homogeneous group.

| Variable                      | common (n=9116)    | outliers (n=480)   | P-value      | Adjusted P-value |
|:------------------------------|-------------------:|-------------------:|-------------:|-----------------:|
| lof                           |      1.06 (0.06)   |       1.49 (0.24)  | 1.43e-151    | 4.86e-150        |
| duration                      |     14.64 (9.81)   |      27.76 (20.1)  | 1.15e-38     | 1.96e-37         |
| age                           |     47.15 (13.14)  |      56.51 (19.53) | 5.84e-23     | 6.62e-22         |
| pem_n                         |      4.64 (0.64)   |        4.00 (1.61) | 1.07e-16     | 9.10e-16         |
| fatigue_n                     |      4.53 (0.69)   |        4.14 (1.39) | 1.30e-09     | 8.84e-09         |
| muscleorjointpain_n           |      3.60 (1.17)   |        3.24 (1.48) | 2.34e-07     | 1.33e-06         |
| multiplechemicalsensitivities_n |   13.23 %        |        21.67 %     | 7.68e-07     | 3.73e-06         |
| relatives_n                   |     14.56 %        |        22.08 %     | 1.78e-05     | 7.56e-05         |
| age_ill                       |     32.52 (11.98)  |       28.74 (23.4) | 4.96e-04     | 1.87e-03         |
| interstitialcystitis_n        |      2.97 %        |         6.04 %     | 6.31e-04     | 1.95e-03         |
| stable_n                      |     24.03 %        |        17.29 %     | 6.20e-04     | 1.95e-03         |
| deteriorating_n               |     45.39 %        |        52.71 %     | 1.90e-03     | 5.38e-03         |
| infection_n                   |     50.35 %        |        43.54 %     | 3.68e-03     | 9.62e-03         |
| sleepproblems_n               |      3.85 (1.08)   |        3.70 (1.36) | 1.91e-02     | 4.06e-02         |
| no_trigger_n                  |     14.12 %        |        18.12 %     | 1.90e-02     | 4.06e-02         |
| other_n                       |     14.12 %        |        18.12 %     | 1.90e-02     | 4.06e-02         |

<p align="left">
  <em> Table 6. Comparison between patients with LOF above the 0.95 quantile (outliers) and the others (common). </em>
</p>

### Machine Learning

#### Supervised classification of severity

The performances of the logistic and Naive Bayes models in predicting the severity of the patients from the other variables are documented in Figures 5 and 6, respectively. The logistic classifier shows a better profile: using as a cut-off a predicted probability of 0.75, it reaches a sensitivity and a specificity of 0.75.

<img src="https://github.com/user-attachments/assets/9ef26982-d128-4b67-a0cd-c388ae2131c4" width="800" alt="Supervised_log_Resampling_severity_n">
<p align="left">
  <em>Figure 5. Performances of a logistic regression, with severity as the outcome variable (its value is one for patients who are housebound or worse, zero otherwise). For each value of the cut-off of the predicted probability used to classify the subject, resampling was randomly performed 40 times; then, mean sensitivity and specificity were calculated and plotted. Left: sensitivity, specificity, accuracy, type I error (ERI), and type II error (ERII) are plotted as a function of the predictive-probability cut-off. Right: Sensitivity is plotted as a function of false positive rate. </em>
</p>

<img src="https://github.com/user-attachments/assets/9245da43-510b-4299-80d2-069586fa4ebf" width="800" alt="Supervised_NB_Resampling_severity_n">
<p align="left">
  <em>Figure 6. Performances of a Naive Bayes supervised classification, with severity as exit variable. For details, see the caption of Figure 5. </em>
</p>

The significant coefficients of the logistic regression applied to all the available data are reported in Table 8. We see that the significant positive predictors of severity are fatigue, dizziness, PEM, sensitivity,  sleep problems, and the presence of relatives with the same disease. On the other hand, a course of the disease that is stable, fluctuating or improving negatively regresses with severity. Also, starting the illness after a traumatic life event predicts a less severe profile. The Area Under the Curve of this model is 0.83 (see ROC_log_severity_n.jpg in the file section).

| Variable          | Estimate | SE         | Z       | p-value    |Sign.    |
|:------------------|---------:|-----------:|--------:|-----------:|:--------|
| (Intercept)       |  -5.6934 |     0.3066 | -18.569 | < 2e-16    | ***     |
| fatigue_n         |   0.9764 |     0.0420 |  23.259 | < 2e-16    | ***     |
| dizzness_n        |   0.2307 |     0.0256 |   9.005 | < 2e-16    | ***     |
| sensitivity_n     |   0.2238 |     0.0278 |   8.045 | 8.63e-16   | ***     |
| pem_n             |   0.3224 |     0.0391 |   8.236 | < 2e-16    | ***     |
| stable_n          |  -1.1637 |     0.0750 | -15.527 | < 2e-16    | ***     |
| fluctuating_n     |  -1.1767 |     0.0744 | -15.812 | < 2e-16    | ***     |
| improving_n       |  -2.0894 |     0.1103 | -18.936 | < 2e-16    | ***     |
| trauma_n          |  -0.3639 |     0.1335 |  -2.727 | 6.39e-03   | **      |
| sleepproblems_n   |   0.0663 |     0.0271 |   2.447 | 1.44e-02   | *       |
| relatives_n       |   0.1975 |     0.0825 |   2.394 | 1.67e-02   | *       |

<p align="left">
  <em> Table 8. Significant coefficients of the logistic regression between a binary variable for severity (1 if house-bound or worse, 0 otherwise) and all the available variables. </em>
</p>

#### Supervised classification of outliers

I also trained two models (logistic regression and Naive Bayes) to predict LOF. Their performances in predicting outliers (LOF>0.95th percentile) are reported in Figures 7 and 8, respectively. Naive-Bays outperforms logistic regression in this case.

<img src="https://github.com/user-attachments/assets/69a52aa7-3e4e-4578-9843-b5262ae122b7" width="800" alt="Supervised_log_Resampling_lof">
<p align="left">
  <em>Figure 7. Performances of a logistic regression, with LOF as the outcome variable (its value is one for LOF above the 95th percentile, zero otherwise). For each value of the cut-off of the predicted probability used to classify the subject, resampling was randomly performed 40 times. Next, the mean sensitivity and specificity were calculated and plotted. Left: sensitivity, specificity, accuracy, type I error (ERI), and type II error (ERII) are plotted as a function of the predictive-probability cut-off. Right: Sensitivity is plotted as a function of false positive rate. </em>
</p>

<img src="https://github.com/user-attachments/assets/3a85c51c-e101-4339-9ddd-5dacc029a93d" width="800" alt="Supervised_NB_Resampling_lof">
<p align="left">
  <em>Figure 8. Performances of a Naive Bayes model, with LOF as the outcome variable (its value is one for LOF above the 95th percentile, zero otherwise). See Figure 7 for details. </em>
</p>

The significant coefficients of the logistic regression are indicated below (Table 9). We see that fatigue, sensitivity, pain, and PEM predict non-outlier status. A younger age at first symptoms predicts outlier status, as well as severity. Fluctuations or stability in disease course predict non-outlier status. 

| Variable | Estimate | Std. Error | z value | p-value | Sig |
|---------|---------|-----------|--------|--------|-----|
| fatigue_n | -0.375022 | 0.068515 | -5.474 | 4.41e-08 | *** |
| muscleorjointpain_n | -0.239960 | 0.050479 | -4.754 | 2.00e-06 | *** |
| sensitivity_n | 0.107491 | 0.052813 | 2.035 | 0.04182 | * |
| pem_n | -0.735165 | 0.055192 | -13.320 | <2e-16 | *** |
| interstitialcystitis_n | 0.503695 | 0.232031 | 2.171 | 0.02995 | * |
| severity_n | 0.463162 | 0.079351 | 5.837 | 5.32e-09 | *** |
| age_ill | -0.074715 | 0.004436 | -16.842 | <2e-16 | *** |
| age | 0.093339 | 0.004700 | 19.861 | <2e-16 | *** |
| stable_n | -0.428085 | 0.154072 | -2.778 | 0.00546 | ** |
| fluctuating_n | -0.297612 | 0.140359 | -2.120 | 0.03398 | * |

<p align="left">
  <em> Table 9. Significant coefficients that predict outlier status (1 for outlier, 0 for non-outlier) using a logistic regression. </em>
</p>

#### Supervised classification of deterioration

I trained the logistic and the Naive Bayes models to recognise subjects who deteriorate over time, using as an outcome a binary variable whose value is one when course_n is  5 (fluctuating initially, then mostly deterioration) or 7 (mainly deterioration), zero otherwise. The performances of the two models, evaluated by resampling, are described by Figures 9 and 10.

<img src="https://github.com/user-attachments/assets/12e4d3e3-d56a-4155-ae17-5d0390cf3a9c" width="800" alt="Supervised_log_Resampling_course_n">
<p align="left">
  <em>Figure 9. Performances of a logistic regression, using as an outcome a binary variable whose value is one when course_n is  5 (fluctuating initially, then mostly deterioration) or 7 (mainly deterioration), zero otherwise. For each value of the cut-off of the predicted probability used to classify the subject, resampling was randomly performed 40 times. Next, the mean sensitivity and specificity were calculated and plotted. Left: sensitivity, specificity, accuracy, type I error (ERI), and type II error (ERII) are plotted as a function of the predictive-probability cut-off. Right: Sensitivity is plotted as a function of false positive rate. </em>
</p>

<img src="https://github.com/user-attachments/assets/f796cd9c-842d-48b9-94b1-69c1f721d100" width="800" alt="Supervised_NB_Resampling_course_n">
<p align="left">
  <em>Figure 8. Performances of a Naive Bayes model. See Figure 9 for details. </em>
</p>

The significant coefficients of the logistic regression are in Table 10. Positive predictors of deterioration are pain, fatigue, dizziness, sleep problems, PEM, Fibromyalgia, IBS, asthma, and disease severity. Having relatives with ME is a protective factor against deterioration, as well as a younger age at first symptoms. 

| Variable | Estimate | p-value | Significance |
|---------|---------|--------|-------------|
| Intercept | -9.27 | <2e-16 | *** |
| Relatives with ME | -0.19 | 0.004 | ** |
| Fatigue | 0.29 | 6.62e-14 | *** |
| Muscle or Joint Pain | 0.08 | 0.00021 | *** |
| Dizziness | 0.07 | 0.00032 | *** |
| Sleep Problems | 0.11 | 4.37e-06 | *** |
| Post-Exertional Malaise (PEM) | 0.17 | 7.04e-06 | *** |
| Fibromyalgia | 0.21 | 0.00011 | *** |
| Irritable Bowel Syndrome | 0.17 | 0.00062 | *** |
| Allergies | -0.16 | 0.0016 | ** |
| Asthma | 0.14 | 0.027 | * |
| Disease Severity | 1.03 | <2e-16 | *** |
| Age at Illness Onset | -0.005 | 0.041 | * |
| Age | 0.009 | 0.00020 | *** |
| LOF | 1.37 | 2.84e-09 | *** |

<p align="left">
  <em> Table 10. Significant coefficients that predict deterioration status (1 for deterioration, 0 otherwise) using a logistic regression. </em>
</p>

#### Prediction for three patients

I applied three of the classifiers discussed above to the supplementary patients included in folder [Supplementary Cases](Supplementary%20Cases/) (Table 11). The cut-off for the three classifiers employed is as follows: 75% for the housebound status predictor (logistic, Figure 5), 5% for the outlier classifier (NB model of Figure 8), and 50% for the classifier of deterioration (logistic regression of Figure 9) 

| Patient | Outlier | Prediction | Housebound | Prediction | Deteriorating | Prediction |
|:--------|--------:|-----------:|-----------:|-----------:|--------------:|-----------:|
| PT1     | Yes     | Yes (10%)  | Yes        | Yes (87%)  | No            | No (45%)   |
| PT2     | No      | No (4%)    | Yes        | Yes (83%)  | Yes           | Yes (65%)  |
| PT2     | No      | No (1%)    | Yes        | No (40%)   | No            | Yes (58%)  |

<p align="left">
  <em>Table 2. Predictions about outlier status, deterioration, and severity for the three patients included in folder [Supplementary Cases](Supplementary%20Cases/). </em>
</p>

#### Unsupervised clustering

Partitioning (clustering) of the data into k clusters around medoids, for k from 2 to 7, failed at finding separated clusters (see folder [Unsupervised](Unsupervised/) for details). Principal Component Analysis provides little contribution to this analysis because 80% of the total variance is distributed across 16 principal components (Figure 3).

# Discussion

Monogenic forms of common idiopathic diseases have been documented in Alzheimer's Disease (AD) ([Barber RC 2012](https://pubmed.ncbi.nlm.nih.gov/24278680/)), Parkinson's Disease (PD) ([Girija MS et Kishore A 2025](https://journals.lww.com/aomd/fulltext/2025/01000/clinical_approach_to_monogenic_parkinson_s_disease.1.aspx)), and Frontotemporal Dementia (FTD) ([Barbier M et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30168435/)). Using PD as a paradigmatic case, we observe that despite a substantial overlap in clinical presentation, monogenic forms often display atypical symptoms (due to the specific gene affected) and an early onset. Also, the study of monogenic cases has offered important insights into the mechanisms of the idiopathic form ([Girija MS et Kishore A 2025](https://journals.lww.com/aomd/fulltext/2025/01000/clinical_approach_to_monogenic_parkinson_s_disease.1.aspx)). 

Myalgic Encephalomyelitis/Chronic Fatigue Syndrome (ME/CFS) is a complex disease defined by persistent fatigue, exacerbation of post-exertional symptoms, cognitive impairment and orthostatic intolerance ([Chu L et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30805319/)), with a prevalence of 0.009 ([Lim EJ et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32093722/)). A recent GWAS meta-analysis on 21,560 ME/CFS patients inclunding data from DecodeME, UK Biobank, and Million Veteran Program, revealed a possible involvement of the glutamatergic synapses ([Maccallini 2025](https://github.com/paolomaccallini-hub/MetaME)) ([Maccallini 2026](https://github.com/paolomaccallini-hub/MetaME-GeNet)). On the other hand, only a few attempts have been published to document Mendelian ME/CFS (MME) by next-generation sequencing (NGS), tagging genes NOS3 ([McGarrity S. et al. 2024](https://pubmed.ncbi.nlm.nih.gov/39328057/)) and AKRIC1/AKRIC2  ([Oakley J et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37978513/)). With the planned whole-genome sequencing (WGS) study of more than 18k patients announced by the DecodeME group, this landscape is likely to change drastically, with a sharp increase in well-documented MME cases.

Soon, clinicians will face the challenge of selecting subjects with ME/CFS to send for WGS, along with their family members. We do not yet know what peculiar features to look for in this population. Here, I propose an unsupervised selector based on the Local Outlier Factor (LOF) algorithm ([Breunig MM, 2000](https://dl.acm.org/doi/10.1145/335191.335388)) applied to the raw data from the EMEA Pan-European ME Patient Survey ([Angelsen A et Schei T, 2024](https://www.europeanmealliance.org/emea-pan-european-survey-UK.shtml)).
