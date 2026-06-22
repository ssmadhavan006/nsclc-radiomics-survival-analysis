# Workspace Rules: Quantitative Radiomic Characterization of NSCLC Tumors

This document defines the core operational boundaries, research guidelines, and coding standards for this project. Every decision, script execution, and analysis method must adhere to these rules.

---

## 1. Research Integrity Rules

### Rule 1
Never fabricate results, metrics, findings, significance values, feature counts, or conclusions.

### Rule 2
Every scientific claim must be traceable to:
* extracted data,
* statistical output,
* literature citation,
* or explicit computation.

### Rule 3
Never claim clinical significance when only statistical significance has been established.

### Rule 4
Clearly separate:
* exploratory analysis,
* statistical association,
* predictive modeling,
* causal interpretation.

---

## 2. Dataset Handling Rules

### Rule 5
Treat the TCIA dataset as the single source of truth. Never:
* hardcode labels,
* synthesize patient records,
* invent missing metadata.

### Rule 6
Every patient included in analysis must have:
* valid CT volume,
* valid segmentation,
* successful feature extraction.

Otherwise exclude and log.

### Rule 7
Every excluded patient must appear in `failed_cases.csv` with a clear reason.

---

## 3. DICOM Safety Rules

### Rule 8
Never assume DICOM slices are correctly ordered. Always sort using DICOM metadata (e.g., Image Position Patient Z).

### Rule 9
Never assume segmentation aligns with CT. Always verify:
* spacing,
* origin,
* direction,
* dimensions,
before extraction.

### Rule 10
Abort feature extraction if mask overlap validation fails. Never continue with potentially misaligned masks.

---

## 4. Radiomics Rules

### Rule 11
Radiomic features must only be computed inside the tumor ROI. No whole-image extraction.

### Rule 12
All preprocessing parameters must be explicitly stored. Example:
* `voxel_spacing: [1, 1, 1]`
* `bin_width: 25`
* `interpolator: sitkBSpline`

### Rule 13
Feature extraction settings must be version controlled. No hidden defaults.

### Rule 14
Raw features and processed features must be stored separately. Example:
* `raw_features.csv`
* `processed_features.csv`

---

## 5. Reproducibility Rules

### Rule 15
Every experiment must be reproducible from a single command. Example:
* `python run_pipeline.py`

### Rule 16
All randomness must use fixed seeds. Example:
* `SEED = 42`

### Rule 17
Every generated artifact must include metadata:
* timestamp,
* git commit,
* configuration version.

---

## 6. Statistical Analysis Rules

### Rule 18
No feature-level hypothesis testing without multiple-testing correction. Use:
* Benjamini-Hochberg FDR,
* Bonferroni,
where appropriate.

### Rule 19
Report effect sizes alongside p-values.

### Rule 20
Never rank features solely by p-value.

### Rule 21
Statistical tests must validate assumptions before execution (e.g., normality, variance homogeneity, sample size adequacy).

---

## 7. Survival Analysis Rules

### Rule 22
Use survival time and event status exactly as defined. Never redefine endpoints.

### Rule 23
No patient leakage between training and testing.

### Rule 24
Use concordance index for survival evaluation. Accuracy is not acceptable for survival prediction.

### Rule 25
All Kaplan-Meier curves must include:
* risk table,
* confidence intervals.

---

## 8. Machine Learning Rules

### Rule 26
Classical ML only.
* **Allowed**: Cox PH, Random Survival Forest, Logistic Regression, Random Forest, XGBoost.
* **Forbidden** (unless explicitly approved): CNNs, ResNet, Vision Transformers, Deep Learning.

### Rule 27
Feature selection must occur inside cross-validation. Never before train/test split.

### Rule 28
Prevent leakage at all stages. This is critical.

---

## 9. Software Engineering Rules

### Rule 29
No file may exceed 500 lines unless justified.

### Rule 30
Every function must:
* have type hints,
* docstrings,
* error handling.

### Rule 31
No duplicated logic. Refactor common operations.

### Rule 32
No notebook-only code. Everything must be runnable as scripts/modules.

---

## 10. Logging Rules

### Rule 33
Every pipeline stage must log:
* start,
* finish,
* duration,
* status,
* errors.

### Rule 34
Never suppress exceptions silently.

### Rule 35
Every warning must be actionable.

---

## 11. Visualization Rules

### Rule 36
Every figure must be publication quality. Include:
* title,
* axis labels,
* legends,
* units.

### Rule 37
Never use screenshots as final outputs. Generate figures programmatically.

---

## 12. Project Scope Protection Rules

### Rule 38
The project objective is: **Quantitative Radiomic Characterization of NSCLC Tumors and Their Association with Clinical Outcomes**. Do not drift into:
* cancer detection,
* CT classification,
* medical chatbot,
* dashboard-first development,
* generic AI applications.

### Rule 39
Whenever uncertainty exists, choose:
`Scientific validity` > `Engineering convenience` > `Visual appeal`.

### Rule 40
Before implementing any major component, verify: Does this directly support radiomic feature extraction, tumor characterization, or clinical association analysis? If not, reject it.

---

## 13. Commands and Scripts Execution Rule

### Rule 41 (Important Execution Limit)
If the assistant wants to run any scripts or commands that are expected to take a long time (such as installing packages or training models), the assistant **must** provide the exact command to the user so the user can run it personally and return the output. For small, short-duration run commands (like fast folder listings, file reads, or small validations), the assistant may run them directly on its own.
