band_arr=("u" "g" "r" "i" "z" "y")
len_arr=(5 6 13 15 13 16)
band_index=3

#----------------------------
band="${band_arr[$band_index]}"
visit_max="${len_arr[$band_index]}"

#============================
sed "s/BAND_INPUT/'${band}'/g" lib/tools_template.py > lib/tools.py

jupyter nbconvert --execute --to notebook --inplace step0*.ipynb


for visit_index in {1..${visit_max}}; do

    echo ""
    echo "========== visit index: ${visit_index} =========="
    
    jupyter nbconvert --execute --to notebook --inplace step1*.ipynb
    jupyter nbconvert --execute --to notebook --inplace step2*.ipynb

    bash cp.sh ${band} ${visit_index}

done

#----------------------------

#mv stamp stamp_${band}
