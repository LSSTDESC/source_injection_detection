band="$1"
visit_index="$2"
field="ecdfs"
version="v2.0.2"

tag="${version}_${field}_${band}_${visit_index}"


#------------------------------
folder="${HOME}/shared_lagn_injection/${tag}/"

[ ! -d $folder ] && mkdir $folder 

#------------------------------
mv catalog/*tab*.csv $folder

mv catalog/inj_catalog*visit*fits $folder
# No need to use it
[ "${visit}" == "1" ] && mv catalog/inj_catalog*template*fits $folder

#------------------------------
mv fig/triple_*png $folder
mv fig/dual_*png $folder

mv cutout $folder

#------------------------------
#rm fig/template*fits 

rm fig/visit*fits 
rm fig/diff*fits

rm fig/injected_visit*fits 
#rm fig/injected_temp*fits 
rm fig/injected_diff*fits
rm fig/kernel*


#------------------------------
rm -r catalog
