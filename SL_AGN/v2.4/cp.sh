band="$1"
visit_index="$2"
field="ecdfs"
version="v1.3"

tag="${version}_${field}_${band}_${visit_index}"


#------------------------------
folder="${HOME}/shared_lagn_injection/${tag}/"

[ ! -d $folder ] && mkdir $folder 

#------------------------------
mv catalog/*tab*.csv $folder

mv catalog/inj_catalog*visit*fits $folder
mv catalog/inj_catalog*template*fits $folder

#------------------------------
mv fig/triple_*png $folder
mv fig/dual_*png $folder
#mv fig/injected_*fits $folder

mv cutout $folder

#------------------------------
#rm fig/template*fits 
#rm fig/visit*fits 
#rm fig/diff*fits

#rm fig/injected_visit*fits 
#rm fig/injected_diff*fits

#------------------------------
#mv catalog catalog_${tag}
#mv fig fig_${tag}
rm -r catalog fig
