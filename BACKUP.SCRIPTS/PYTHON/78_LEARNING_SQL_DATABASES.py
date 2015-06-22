#LEARNING SQL DATABASES
import MySQLdb as MQ
from FUNCTIONS import unique

#Create connection object that represents database
conn_obj=MQ.connect(host="localhost",user="root", passwd="mysql", db="TESTING")

cursor=conn_obj.cursor()

cursor.execute("SELECT NAME FROM CSA_catalyticCofactors")

print cursor.fetchall()

cursor.execute("SHOW TABLES FROM TESTING")
for table in cursor.fetchall():
    print list(table)
    
cursor.execute("SHOW FIELDS FROM CSA_catalyticCofactors")
print cursor.fetchall()

cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_name='CSA_catalyticCofactors'")
print cursor.fetchall()

cursor.execute("SELECT PDBID, CHAIN,NUMBER, UPNUMBER , TYPE FROM CSA_CATALYTICRESIDUES")
RESULTS= [list(X) for X in cursor.fetchall()]
print len(RESULTS)
RESULTS=unique(RESULTS)
print RESULTS
print len(RESULTS)

RESULTS_NULL=[X for X in RESULTS if X[3]==None]
print len(RESULTS_NULL), RESULTS_NULL

