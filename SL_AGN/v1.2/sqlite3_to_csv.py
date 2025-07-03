# Convert a sqlite3 database file to a csv table
# We only select one table from the sqlite3 file
# Note we can also do it through SQL command line
# Section 7.6. Export to CSV in https://www.sqlite.org/cli.html 
# https://pipelines.lsst.io/getting-started/dc2-guide.html



#=======================
from astropy.table import Table
import sys
import sqlite3
import pandas as pd



#=======================
if len(sys.argv) != 4:
    print("Usage: python this.py in.sqlite3 table_name out.csv")
    print("Exiting...")
    sys.exit(1)

in_filename = sys.argv[1]
table_name = sys.argv[2]
out_filename = sys.argv[3]



#-----------------------
in_sqlite3 = sqlite3.connect(in_filename)

# Print a list of all tables
cursor = in_sqlite3.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
print(tables)
cursor.close()

# Get a table
table_sqlite3 = pd.read_sql_query("select * from %s;"%table_name, in_sqlite3)
in_sqlite3.close()

# Save the table
#table_astropy = Table.from_pandas(table_sqlite3)
#table_astropy.write(out_filename, format="ascii.csv", overwrite=True)
table_sqlite3.to_csv(out_filename)
