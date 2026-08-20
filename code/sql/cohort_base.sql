  SELECT
    i.stay_id, i.hadm_id, i.subject_id, i.intime, i.outtime, i.los,
    p.gender, p.anchor_age,
    a.admission_type, a.insurance, a.hospital_expire_flag,
    MAX(CASE WHEN (d.icd_version = 9  AND d.icd_code LIKE '428%')
             OR (d.icd_version = 10 AND d.icd_code LIKE 'I50%') THEN 1 ELSE 0 END) AS chf,
    MAX(CASE WHEN (d.icd_version = 9  AND d.icd_code LIKE '250%')
             OR (d.icd_version = 10 AND (d.icd_code LIKE 'E10%' OR d.icd_code LIKE 'E11%'
                                       OR d.icd_code LIKE 'E12%' OR d.icd_code LIKE 'E13%'
                                       OR d.icd_code LIKE 'E14%')) THEN 1 ELSE 0 END) AS diabetes,
    MAX(CASE WHEN (d.icd_version = 9  AND d.icd_code IN ('491','492','494','496'))
             OR (d.icd_version = 10 AND d.icd_code LIKE 'J44%') THEN 1 ELSE 0 END) AS copd,
    MAX(CASE WHEN (d.icd_version = 9  AND (d.icd_code LIKE '584%' OR d.icd_code LIKE '585%' OR d.icd_code LIKE '586%'))
             OR (d.icd_version = 10 AND (d.icd_code LIKE 'N17%' OR d.icd_code LIKE 'N18%'
                                       OR d.icd_code LIKE 'N19%')) THEN 1 ELSE 0 END) AS renal_failure,
    MAX(CASE WHEN (d.icd_version = 9  AND d.icd_code LIKE '401%')
             OR (d.icd_version = 10 AND (d.icd_code LIKE 'I10%' OR d.icd_code LIKE 'I11%'
                                       OR d.icd_code LIKE 'I12%' OR d.icd_code LIKE 'I13%'
                                       OR d.icd_code LIKE 'I15%')) THEN 1 ELSE 0 END) AS hypertension
  FROM icustays i
  LEFT JOIN patients p      ON i.subject_id = p.subject_id
  LEFT JOIN admissions a    ON i.hadm_id = a.hadm_id
  LEFT JOIN diagnoses_icd d ON i.hadm_id = d.hadm_id
  GROUP BY i.stay_id, i.hadm_id, i.subject_id, i.intime, i.outtime, i.los,
           p.gender, p.anchor_age, a.admission_type, a.insurance, a.hospital_expire_flag
  ORDER BY i.stay_id
