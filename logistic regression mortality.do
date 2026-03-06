/* ============================================================
   PROJECT:     Logistic Regression Analysis of 3-Year Mortality
   DATASET:     mortality.dta
   COURSE:      Statistical Methods in Epidemiology
                London School of Hygiene & Tropical Medicine (LSHTM)
   AUTHOR:      Uchechi Lever
   DATE:        2025

   DESCRIPTION:
   This do file analyses predictors of 3-year all-cause mortality
   in a cohort study using logistic regression. The study examines
   a population-based cohort followed for a fixed period of 3 years.

   KEY EXPOSURES:
     vimp    - Visual impairment (binary: 1=impaired, 0=unimpaired)
     mfgrp   - Microfilarial load grouped into 4 levels
                 0 = uninfected (reference)
                 1 = <10 mf/mg
                 2 = 10-49 mf/mg
                 3 = >=50 mf/mg

   OUTCOME:
     died    - 3-year mortality (binary: 1=died, 0=alive)

   COVARIATES:
     agegrp  - Age group (categorical)
     sex     - Sex (binary)
     bmi     - Body mass index (continuous)
     education - Education level (categorical)

   METHODS:
   - Descriptive statistics and tabulations at baseline
   - Crude odds ratios via tabodds and mhodds
   - Logistic regression: crude and adjusted models
   - Indicator (dummy) variable approach for categorical exposures
   - Likelihood ratio tests for categorical variables
   - Confounding assessment (vimp adjusted for age)

   OUTPUT:
   - Baseline characteristics table
   - Odds ratios with 95% CIs for vimp and mfgrp
   - Likelihood ratio test statistics
   - Adjusted model comparing crude vs adjusted OR for vimp

   NOTE:
   This analysis treats the data as a cohort study with fixed
   follow-up time and uses odds as an approximation of risk.
   All analyses conducted in Stata 17.
   ============================================================ */


/* ============================================================
   SECTION 0: HOUSEKEEPING
   ============================================================ */

   clear all
   set more off
   capture log close

   * Set working directory — update this path to your local folder
   * cd "your/path/here"

   * Open log file to record all output
   log using "logistic_regression_mortality_output.log", replace text

   * Load dataset
   use mortality, clear


/* ============================================================
   SECTION 1: DATA EXPLORATION
   Understand the structure, variable types, and missing data
   before any analysis. This is essential good practice.
   ============================================================ */

   * Overview of dataset structure
   describe

   * Summary statistics for all variables
   summarize

   * Inspect first 20 rows to sense-check data
   list in 1/20

   * Check outcome variable coding
   tab died, missing

   * Check key exposure variables for missing data
   tab vimp, missing
   tab mfgrp, missing     /* Note: mfgrp has 93 missing values */
   tab agegrp, missing


/* ============================================================
   SECTION 2: BASELINE CHARACTERISTICS
   Describe the cohort at enrolment (before follow-up).
   This forms the basis of the "Table 1" in a results section.
   ============================================================ */

   * Age distribution
   summarize age, detail

   * Sex distribution
   tab sex

   * Visual impairment prevalence at baseline
   tab vimp

   * Microfilarial infection by group
   tab mfgrp

   * Age group distribution
   tab agegrp

   * BMI summary
   summarize bmi, detail

   * Education level
   tab education

   * Cross-tabulation: visual impairment by sex
   tab vimp sex, row


/* ============================================================
   SECTION 3: CRUDE ANALYSIS — VISUAL IMPAIRMENT AND DEATH
   Examine the unadjusted association between visual impairment
   (vimp) and 3-year mortality using tabodds and mhodds.
   ============================================================ */

   * Step 3a: Simple cross-tabulation of died by vimp
   * Row percentages show risk of death within each vimp group
   tab died vimp, row col

   * Step 3b: Tabodds — displays odds of death at each level of vimp
   * Useful first step before formal regression
   tabodds died vimp

   * Step 3c: Mhodds — Mantel-Haenszel odds ratio for vimp vs died
   * Gives a single crude OR comparing visually impaired vs unimpaired
   mhodds died vimp


/* ============================================================
   SECTION 4: LOGISTIC REGRESSION — VISUAL IMPAIRMENT AND DEATH
   Formal logistic regression models for the association between
   visual impairment and 3-year mortality.
   ============================================================ */

   * Step 4a: Logit command — output on the LOG ODDS scale
   * Useful for understanding iteration process and log likelihood
   logit died vimp

   /*
   Key output to note:
   - Log likelihood at convergence
   - LR chi2 statistic: tests H0 that vimp is not associated with death
   - Number of observations
   */

   * Step 4b: Logistic command — output on the ODDS RATIO scale
   * This is the preferred format for reporting results
   logistic died vimp

   /*
   Interpretation guide:
   - Coefficient on vimp = OR for death comparing visually impaired
     to visually unimpaired individuals
   - _cons = baseline odds of death (visually unimpaired group)
   - z statistic = Wald test for each parameter
   - 95% CI derived from standard error of the LOG odds ratio
   */

   * Step 4c: Equivalent syntax using the 'or' option
   logit died vimp, or


/* ============================================================
   SECTION 5: CATEGORICAL EXPOSURE — MICROFILARIAL LOAD (mfgrp)
   mfgrp has 4 ordered categories requiring indicator variables.
   The reference group is mfgrp=0 (uninfected).
   ============================================================ */

   * Step 5a: Cross-tabulation — column percentages appropriate here
   * because exposure (mfgrp) determines column grouping
   tab mfgrp died, col

   * Step 5b: Logit model with indicator variables (log odds scale)
   * The i. prefix automatically generates indicator (dummy) variables
   logit died i.mfgrp

   * Step 5c: Logistic model with indicator variables (OR scale)
   * 'base' option shows the reference category explicitly in output
   logistic died i.mfgrp, base

   /*
   Expected odds ratios (vs uninfected reference group mfgrp=0):
     mfgrp=1 (<10 mf/mg):    OR = 1.69
     mfgrp=2 (10-49 mf/mg):  OR = 1.46
     mfgrp=3 (>=50 mf/mg):   OR = 2.05

   There are THREE Wald p-values (one per OR).
   There is ONE Likelihood Ratio Test p-value (overall association).
   */

   * Step 5d: Inspect indicator variables created internally by Stata
   * Useful to understand how i. prefix works
   list id mfgrp i.mfgrp in 1/25


/* ============================================================
   SECTION 6: LIKELIHOOD RATIO TEST — MICROFILARIAL LOAD
   Tests the overall null hypothesis that mfgrp is not associated
   with death (simultaneously tests all 3 parameters).

   IMPORTANT: mfgrp has 93 missing values. Both the null and
   alternative models must use the SAME observations to ensure
   a valid comparison of log likelihoods.
   ============================================================ */

   * Step 6a: Check missing data in mfgrp
   tab mfgrp, missing
   /* Confirms 93 missing records */

   * Step 6b: Fit alternative model (with mfgrp)
   * Missing mfgrp records are automatically dropped — n=4205
   logistic died i.mfgrp
   estimates store A          /* Store log likelihood as "A" */

   * Step 6c: Fit null model (without mfgrp) restricted to SAME sample
   * Use 'if mfgrp!=.' to exclude the 93 missing mfgrp records
   * This ensures both models use n=4205 observations
   logistic died if mfgrp!=.
   estimates store B          /* Store log likelihood as "B" */

   * Step 6d: Likelihood ratio test comparing models A and B
   lrtest B A

   /*
   Expected output:
     LR chi2(3) = 6.70
     Prob > chi2 = 0.0822

   Interpretation: There is weak evidence (p=0.08) of an overall
   association between microfilarial load and 3-year mortality.
   */


/* ============================================================
   SECTION 7: CONFOUNDING ASSESSMENT — AGE AS A CONFOUNDER
   Age has been identified as a potential confounder of the
   relationship between visual impairment and death (causal graph).
   We assess confounding by comparing crude and adjusted ORs.
   ============================================================ */

   * Step 7a: Crude model — vimp only (repeated from Section 4)
   logistic died vimp

   * Step 7b: Age alone — assess age-death relationship
   * Uses indicator variables for categorical agegrp
   logistic died i.agegrp, base

   /*
   Interpretation: This tells us whether age is independently
   associated with death — a prerequisite for it to be a confounder.
   */

   * Step 7c: Adjusted model — vimp + agegrp simultaneously
   * If OR for vimp changes materially, age is a confounder
   logistic died i.vimp i.agegrp, base

   /*
   Compare the OR for vimp from:
     - Crude model  (Section 4):   unadjusted OR
     - Adjusted model (this step): OR after controlling for age

   A meaningful change in the vimp OR suggests age confounds
   the vimp-death association.
   */


/* ============================================================
   SECTION 8: RESULTS SUMMARY
   Consolidated display of key model results for reporting.
   ============================================================ */

   * Final crude model: visual impairment
   di _newline "--- CRUDE MODEL: Visual Impairment ---"
   logistic died vimp

   * Final crude model: microfilarial load
   di _newline "--- CRUDE MODEL: Microfilarial Load ---"
   logistic died i.mfgrp, base

   * Final adjusted model: vimp adjusted for age
   di _newline "--- ADJUSTED MODEL: vimp + agegrp ---"
   logistic died i.vimp i.agegrp, base

   * Close log
   log close


/* ============================================================
   END OF DO FILE

   INTERPRETATION NOTES FOR GITHUB READERS:
   ------------------------------------------
   This analysis is based on a cohort study examining 3-year
   all-cause mortality in relation to two key exposures:

   1. VISUAL IMPAIRMENT (vimp):
      Visually impaired individuals had higher odds of dying
      within 3 years compared to the unimpaired (crude OR ~1.9,
      p<0.001). This association persisted after adjustment
      for age, suggesting visual impairment is an independent
      predictor of mortality.

   2. MICROFILARIAL LOAD (mfgrp):
      Higher microfilarial load was associated with increased
      odds of death, with the highest load group (>=50 mf/mg)
      showing OR ~2.05 vs uninfected. However, the overall
      LR test (p=0.08) suggests weak evidence at conventional
      thresholds, warranting cautious interpretation.

   3. CONFOUNDING BY AGE:
      Age was a strong independent predictor of death. Comparing
      crude and age-adjusted ORs for vimp assesses whether the
      vimp-mortality association is confounded by age.

   SOFTWARE: Stata 17
   DATASET:  mortality.dta (LSHTM SME course dataset)
             Note: raw data not included in this repository.
             A simulated version with equivalent structure can
             be generated using the seed script: simulate_data.do
   ============================================================ */
