setup lsst_distrib

#============================
[ ! -d "catalog" ] && mkdir "catalog"
[ ! -d "fig" ] && mkdir "fig"


python set_stamp.py

#tar -cvzf fig.tgz fig
