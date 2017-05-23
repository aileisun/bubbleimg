SELECT
  main.object_id, main.ra, main.dec, main.patch_id, main.tract, main.patch, main.patch_s, main.parent_id, main.deblend_nchild, main.detect_is_patch_inner, main.detect_is_tract_inner, main.detect_is_primary
FROM
  {layer}.forced AS main
  
WHERE
  (main.detect_is_tract_inner = 't' AND main.detect_is_patch_inner = 't' AND main.detect_is_primary ='t') AND (coneSearch(coord, {ra}, {dec}, {radius}))

LIMIT
  10