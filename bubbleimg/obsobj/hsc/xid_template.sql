SELECT
  object_id, 
  ra, 
  dec, 
  patch_id, 
  tract, 
  patch, 
  patch_s, 
  parent_id, 
  deblend_nchild, 
  detect_is_patch_inner, 
  detect_is_tract_inner, 
  detect_is_primary

FROM
  {rerun}.forced 
  
WHERE
  (detect_is_tract_inner = 't' 
  	AND detect_is_patch_inner = 't' 
  	AND detect_is_primary ='t') 
  AND (coneSearch(coord, {ra}, {dec}, {radius}))

LIMIT
  10