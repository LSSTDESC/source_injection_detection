band_arr=("u" "g" "r" "i" "z" "y")
len_arr=(5 6 13 15 13 16)
band_index=0

#----------------------------
band="${band_arr[$band_index]}"
iter_max="${len_arr[$band_index]}"
echo "band: ${band}"
echo "iter_max: ${iter_max}"

#============================
sed "s/BAND_INPUT/'${band}'/g" lib/tools_template.py > lib/tools.py

jupyter nbconvert --execute --to notebook --inplace step0*.ipynb


#for iter_index in {1..2}; do
#for iter_index in $(seq 3 ${iter_max}); do
for iter_index in $(seq 1 ${iter_max}); do

    echo ""
    echo "========== iter index: ${iter_index} =========="

    jupyter nbconvert --execute --to notebook --inplace step1*.ipynb

    # Append record with iter_index
    jq -c ". + {\"iter_index\": ${iter_index}}" numbers.json >> numbers_log.jsonl

    jupyter nbconvert --execute --to notebook --inplace step2*.ipynb

    bash cp.sh ${band} ${iter_index}

done

#----------------------------
rm fig/*temp*fits 
#mv stamp stamp_${band}
