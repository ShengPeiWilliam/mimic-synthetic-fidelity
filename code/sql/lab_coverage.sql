  WITH labs AS (
    SELECT
      itemid,
      CASE
        WHEN LOWER(label) LIKE '%creatinine%'        THEN 'Creatinine'
        WHEN LOWER(label) LIKE '%glucose%'           THEN 'Glucose'
        WHEN LOWER(label) LIKE '%hemoglobin%'        THEN 'Hemoglobin'
        WHEN LOWER(label) LIKE '%white blood cell%'  THEN 'WBC'
        -- exclude lactate dehydrogenase, a different assay
        WHEN LOWER(label) LIKE '%lactate%'
             AND LOWER(label) NOT LIKE '%dehydrogenase%' THEN 'Lactate'
        ELSE NULL
      END AS lab
    FROM d_labitems
    WHERE fluid = 'Blood'
  ),
  total AS (
    SELECT COUNT(*) AS n_stays FROM icustays
  ),
  covered AS (
    SELECT l.lab, COUNT(DISTINCT i.stay_id) AS n_with
    FROM icustays i
    JOIN labevents e ON i.hadm_id = e.hadm_id
    JOIN labs l      ON e.itemid = l.itemid
    WHERE l.lab IS NOT NULL
      AND e.charttime BETWEEN i.intime AND i.intime + INTERVAL 24 HOUR
      AND e.valuenum IS NOT NULL
    GROUP BY l.lab
  )
  SELECT
    c.lab,
    c.n_with,
    t.n_stays,
    ROUND(100.0 * c.n_with / t.n_stays, 1) AS coverage_pct
  FROM covered c, total t
  ORDER BY coverage_pct DESC
