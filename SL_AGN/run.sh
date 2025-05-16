# Running on the Rubin Science Platform

#============================
setup lsst_distrib


#============================
[ ! -d "catalog" ] && mkdir -p "catalog"
[ ! -d "fig" ] && mkdir -p "fig"

#python set_stamp.py

#python inject_source.py

#python diffim_detect.py

#python check_completeness.py

#apdb-cli create-sql "sqlite:///apdb.db" apdb_config.py
python associate.py

#============================
#python sqlite3_to_csv.py apdb.db DiaObject DiaObject.csv
#python sqlite3_to_csv.py apdb.db DiaSource DiaSource.csv
#python sqlite3_to_csv.py apdb.db SSObject SSObject.csv
#python sqlite3_to_csv.py apdb.db metadata metadata.csv
#python sqlite3_to_csv.py apdb.db DiaForcedSource DiaForcedSource.csv
#python sqlite3_to_csv.py apdb.db DiaObject_To_Object_Match DiaObject_To_Object_Match.csv
#rename.ul .csv _empty.csv *.csv

#============================
#tar -cvzf fig.tgz fig
