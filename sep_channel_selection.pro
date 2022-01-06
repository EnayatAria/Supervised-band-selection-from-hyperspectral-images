; 
pro SEP_channel_selection
  ;First restore all the base save files.
  envi, /restore_base_save_files
  
  ;Initialize ENVI and send all errors and warnings to the file batch.txt
  envi_batch_init, log_file='batch.txt'
  
  ;Open the input file after removing the noisy channels
  envi_open_file, $
  ;  'D:\shosseiniaria\Enayat documents\PhD Thesis\chapters\separability\data\Salinas\190_bands',$
  'D:\shosseiniaria\Enayat documents\PhD Thesis\chapters\separability\data\IndPin_178b_Qc_P',$
    r_fid=fid
  if (fid eq -1) then begin
    envi_batch_exit
    return
  endif
  
  ;Get the dimensions and # bands for the input file.
  envi_file_query, fid, dims=dims, nb=nb, wl=wl
  
  ;Set the POS keyword to calculate statistics for all data(spectrally) in the file.
  pos=lindgen(nb)
  
  
  
  ;restore Region of Interest- the place of two given classes
  envi_restore_rois, $
  ;  'D:\shosseiniaria\Enayat documents\PhD Thesis\chapters\separability\data\Salinas\16_classes.roi'
   'D:\shosseiniaria\Enayat documents\PhD Thesis\chapters\separability\data\13classes.roi'
  roi_ids = envi_get_roi_ids(fid=fid, $
    roi_colors=roi_colors, roi_names=class_names)
  class_names = ['Unclassified', class_names]
  num_classes = n_elements(roi_ids)
  npts_mat=intarr(num_classes)
  max_pix=0
  for i=0, num_classes-1 do begin
    ENVI_GET_ROI_INFORMATION, roi_ids(i), npts=npts
    npts_mat(i)=npts
    if npts gt max_pix then max_pix=npts
  endfor
  data=intarr(nb,max_pix,num_classes)
  for i=0, num_classes-1 do begin
    a=envi_get_roi_data(roi_ids[i], fid=fid, pos=pos)
    neg_pix=where (a le 0,c)
    if c gt 0 then a(neg_pix)=1
    ENVI_GET_ROI_INFORMATION, roi_ids(i), npts=npts
    data(*,0:npts-1,i)=temporary(a)
  endfor
  data=transpose(temporary(data))
  ;data=data(*,*,20:39)
  ;nb=20 
  ; computing variance for a band
;  b26=data(*,*,31)
;  ;data1=mean(data(*,*,29:35),dimension=3)
; 
;  a1=0
;  for i=0, num_classes-1 do begin
;    b1=variance(b26(i,0:npts_mat(i)-1))
;   ; b2=variance(data2(i,0:npts_mat(i)-1))
; 
;    print,i+1, b1 ;, b2
;    a1=a1+b1 ; &  a2=a2+b2  
;  endfor
;   print, a1 ;,  a2;,  a3 , a4
;   ;print, mean([a1, a2, a3, a4],/double)
     
  ;External file
  openw, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=100
  printf,lun, 'No.bands  location  wavelength   Divergence '
  free_lun, lun
  
  openw, lun, 'CH_selec_SFFS_elemination.txt',/get_lun,width=100
  printf,lun, 'No.bands  location  Divergence '
  free_lun, lun
  ; npts_mat1=rebin(npts_mat, num_classes,2)
   res_struct={final,res:10,ED:0.0D,sep_mat:dblarr(nb-1)}
 result1=dblarr(3,201)
 k=0
  ;**********************************************************************************
  ;test
 ; data=data(1:2,*,*)
 ; num_classes=2
 ; npts_mat=npts_mat(1:2)
 ; nb=nb
 ; ********************************************************************************
 ; SFS only
; WHILE (1) DO BEGIN
;   PRINT, "K=", K
;   
;   
;   ANS=SFS(DATA,NUM_CLASSES,NPTS_MAT,NB,RES_STRUCT,RESULT1)
;   RESULT1(K)=ANS.RES
; 
;   OPENU, LUN, 'CH_SELEC_SEP_RESULT.TXT',/GET_LUN,WIDTH=200, /APPEND
;   PRINTF, LUN, K+1, ANS.RES, WL(ANS.RES-1),ANS.ED
;   FREE_LUN,LUN
;    IF K EQ 30 THEN BREAK
;    K=K+1
;       ENDWHILE
;  **********************************************************************************     
 
  kk=0
  while (1) do begin
    print, "k=", k
   
   
      ans=SFS(data,num_classes,npts_mat,nb,res_struct,result1) 
      kk=kk+1
      result1(0,kk-1)=kk
      result1(1,kk-1)=ans.res
      result1(2,kk-1)=ans.ED
    
  openu, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=200, /append
  printf, lun, kk, ans.res, wl(ans.res-1),ans.ED
  free_lun,lun
  kkk=kk
  if (kkk eq 30) then break
;  ; call SFFS
   if kk ge 3 then begin 
    ans=SFFS(data,num_classes,npts_mat,nb,res_struct,result1)
     a=where(result1(1,*), c)
  if ans.res ne result1(1,a(c-1)) then begin
   b=where(result1(1,*) eq ans.res)
    KK=KK-1
   openu, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=200, /append
   printf, lun, kk, string(result1(1,b),format='(I8)') + '   is removed', ans.ED
   free_lun,lun
   openu, lun, 'CH_selec_SFFS_elemination.txt',/get_lun,width=100,/append
   printf,lun, result1(*,0:c-1)
    printf,lun,'**********************************************'
   free_lun, lun
   
   result1(*,b)=0
  result1=update_sep_ch_select(data,num_classes,npts_mat,nb,res_struct,result1,b)
  openu, lun, 'CH_selec_SFFS_elemination.txt',/get_lun,width=100,/append
  printf,lun, result1(*,0:c-1)
  free_lun, lun
   
   ; Step 3
   while 1 do begin
       if kk eq 2 then break
    ans=SFFS(data,num_classes,npts_mat,nb,res_struct,result1)
    a=where(result1(1,*), c)
     if ans.res ne result1(1,a(c-1)) then begin
       if ans.ED gt result1(2,c-2) then begin
        b=where(result1(1,*) eq ans.res)
          KK=KK-1
          openu, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=200, /append
          printf, lun, kk, result1(1,b), "is removed", ans.ED
          free_lun,lun
          
          openu, lun, 'CH_selec_SFFS_elemination.txt',/get_lun,width=100,/append
          printf,lun, result1(*,0:c-1)
          printf,lun,'**********************************************'
          free_lun, lun
           result1(*,b)=0
           result1=update_sep_ch_select(data,num_classes,npts_mat,nb,res_struct,result1,b)
           
           openu, lun, 'CH_selec_SFFS_elemination.txt',/get_lun,width=100,/append
           printf,lun, result1(*,0:c-1)
           free_lun, lun

       endif else break
     endif else break
     endwhile
    ;**************************************** 
     
     
  endif 
 endif
 
   if (k eq 200) then break
    if (kk eq 31) then break
  k=k+1
endwhile
;********************************
 a=where(result1, c)
 openu, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=200, /append
printf, lun,  result1(a)
free_lun,lun
; computing the separability of the final set

a=where(result1, c)
arrange=result1(a)-1
data1=data(*,*,arrange)
for i=0, c-1 do begin
 data2=data1(*,*,0:i) 
 ; nbb is the number of band in the subset
nbb=i+1
  npts_mat1=rebin(npts_mat, num_classes,i+1)

  
  
   
    cov=dblarr(nbb,nbb, num_classes)
    sep_mat=dblarr(num_classes, num_classes)
    avg=total(data2,2)/npts_mat1
    
    for j=0, num_classes-1 do begin
      cova=data2(j,0:npts_mat1(j,0)-1,*)
      cova=transpose(temporary(cova))
      if nbb gt 1   then cov[*,*,j]=correlate(cova, /covariance) else cov[*,*,j]=variance(cova)
      
      
    endfor
    TDavg=0 & J_Mavg=0
    
    for t=0, num_classes-1 do begin
      for j=t+1, num_classes-1 do begin
      
      
        Mh=transpose(avg(t,*)-avg(j,*))##invert((cov(*,*,t)+cov(*,*,j))/2)##(avg(t,*)-avg(j,*))
;        ;Compute Battacharyya Distance
        if nbb eq 1 then  B=(Mh/8)+(alog(((cov(*,*,t)+cov(*,*,j))/2)/(sqrt(cov(*,*,t))*sqrt(cov(*,*,j)))))/2 $
        else $
          B=(Mh/8)+(alog((determ((cov(*,*,t)+cov(*,*,j))/2,/double))/(sqrt(determ(cov(*,*,t),/double))*sqrt(determ(cov(*,*,j),/double)))))/2
        ;Compute Divergence
;        if nbb eq 1 then D=(cov(*,*,j)-cov(*,*,t))##(invert(cov(*,*,t))-invert(cov(*,*,j)))/2+$
;          (invert(cov(*,*,t))+invert(cov(*,*,j)))##(avg(j,*)-avg(t,*))##transpose(avg(j,*)-avg(t,*))/2 $
;        else  D=trace((cov(*,*,j)-cov(*,*,t))##(invert(cov(*,*,t))-invert(cov(*,*,j))),/double)/2+$
;          trace((invert(cov(*,*,t))+invert(cov(*,*,j)))##(avg(j,*)-avg(t,*))##transpose(avg(j,*)-avg(t,*)))/2
;        TD=2*(1-exp(-D/8))
        ;   Davg=Davg+D
        J_M= sqrt(2*(1-exp(-B)))
        
       ; sep_mat(t,j)=J_M & sep_mat(j,t)=J_M
        J_Mavg=J_M+J_Mavg
        ;   TDavg=TD+TDavg
      endfor
    endfor
    ; num_classes-1: nth triangular number
    J_Mavg=J_Mavg/(((num_classes-1)^2+num_classes-1)/2)
   
     openu, lun, 'CH_selec_SFFS_sep_result.txt',/get_lun,width=200, /append
    printf, lun,   i+1, arrange(i)+1, wl(arrange(i)),J_Mavg
    free_lun,lun
   
 
endfor

;  envi_batch_exit
end