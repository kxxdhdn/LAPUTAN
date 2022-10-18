function fitsheader_pc2cd,header
  
  h=header

  ;; reply on extast
  extast,header,ast

  sxdelpar,h,'pc1_1'
  sxdelpar,h,'pc1_2'
  sxdelpar,h,'pc2_1'
  sxdelpar,h,'pc2_2'
  sxdelpar,h,'cdelt1'
  sxdelpar,h,'cdelt2'

  putast,h,ast
  
  sxdelpar,h,'pv1_0'
  sxdelpar,h,'pv1_1'
  sxdelpar,h,'pv1_2'
  sxdelpar,h,'pv1_3'
  sxdelpar,h,'pv1_4'

  return,h
end
