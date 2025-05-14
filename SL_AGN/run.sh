setup lsst_distrib


#============================
[ ! -d "catalog" ] && mkdir -p "catalog"
[ ! -d "fig" ] && mkdir -p "fig"

python set_stamp.py

python inject_sourece.py

#python diffim_detect.py

#python associate.py

#============================
#tar -cvzf fig.tgz fig
