import sqlite3

# Create/Connect to database
conn = sqlite3.connect("master.db")
cursor = conn.cursor()

# Read and execute .sql script
with open("database.sql", "r") as f:
    sql_script = f.read()
cursor.executescript(sql_script)

conn.commit()
conn.close()
