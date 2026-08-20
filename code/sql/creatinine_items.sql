  SELECT itemid, label, fluid, category
  FROM d_labitems
  WHERE LOWER(label) LIKE '%creatinine%'
    AND fluid = 'Blood'
  ORDER BY label
