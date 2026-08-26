  SELECT
    i.stay_id,
    MAX(l.valuenum) AS creatinine_max
  FROM icustays i
  JOIN labevents l ON i.hadm_id = l.hadm_id
  WHERE l.itemid IN (50912, 52546, 52024)
    AND l.charttime BETWEEN i.intime AND i.intime + INTERVAL 24 HOUR
    AND l.valuenum IS NOT NULL
  GROUP BY i.stay_id
