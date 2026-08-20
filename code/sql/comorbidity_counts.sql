  SELECT
    CASE
      WHEN icd_version = 9  AND icd_code LIKE '428%' THEN 'CHF'
      WHEN icd_version = 10 AND icd_code LIKE 'I50%'  THEN 'CHF'
      WHEN icd_version = 9  AND icd_code LIKE '250%' THEN 'Diabetes'
      WHEN icd_version = 10 AND (icd_code LIKE 'E10%' OR icd_code LIKE 'E11%'
                               OR icd_code LIKE 'E12%' OR icd_code LIKE 'E13%'
                               OR icd_code LIKE 'E14%') THEN 'Diabetes'
      WHEN icd_version = 9  AND icd_code IN ('491','492','494','496') THEN 'COPD'
      WHEN icd_version = 10 AND icd_code LIKE 'J44%' THEN 'COPD'
      WHEN icd_version = 9  AND icd_code LIKE '584%' THEN 'Renal Failure'
      WHEN icd_version = 9  AND icd_code LIKE '585%' THEN 'Renal Failure'
      WHEN icd_version = 9  AND icd_code LIKE '586%' THEN 'Renal Failure'
      WHEN icd_version = 10 AND (icd_code LIKE 'N17%' OR icd_code LIKE 'N18%'
                               OR icd_code LIKE 'N19%') THEN 'Renal Failure'
      WHEN icd_version = 9  AND icd_code LIKE '401%' THEN 'Hypertension'
      WHEN icd_version = 10 AND (icd_code LIKE 'I10%' OR icd_code LIKE 'I11%'
                               OR icd_code LIKE 'I12%' OR icd_code LIKE 'I13%'
                               OR icd_code LIKE 'I15%') THEN 'Hypertension'
      ELSE NULL
    END AS comorbidity,
    COUNT(DISTINCT d.hadm_id) AS positive_count
  FROM diagnoses_icd d
  WHERE d.hadm_id IN (SELECT DISTINCT hadm_id FROM icustays)
    AND comorbidity IS NOT NULL
  GROUP BY comorbidity
  ORDER BY positive_count
