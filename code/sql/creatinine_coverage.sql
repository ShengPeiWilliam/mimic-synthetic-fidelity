  WITH total AS (
    SELECT COUNT(*) AS n_stays FROM icustays
  ),
  with_creatinine AS (
    SELECT COUNT(DISTINCT i.stay_id) AS n_with_creatinine
    FROM icustays i
    JOIN labevents l ON i.hadm_id = l.hadm_id
    WHERE l.itemid IN (50912, 52546, 52024)
      AND l.charttime BETWEEN i.intime AND i.intime + INTERVAL 24 HOUR
      AND l.valuenum IS NOT NULL
  )
  SELECT
    total.n_stays,
    with_creatinine.n_with_creatinine,
    ROUND(100.0 * with_creatinine.n_with_creatinine / total.n_stays, 1) AS coverage_pct
  FROM total, with_creatinine
