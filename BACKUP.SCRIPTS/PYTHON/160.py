####083114####
#Process Rstudio history_database file

import time

FILE_IN1=open("/Users/jzamalloa/.rstudio-desktop/history_database")
rstudio_history=[[  list(time.gmtime(float(x.split(":",1)[0])/1000)[:4]) , x.split(":",1)[1]] for x in FILE_IN1.read().splitlines()]
FILE_IN1.close()

for i in rstudio_history[:20]:
    print i[0], type(i[0][1]),i[0][0]>2012, i[1]

#Filter for august 2014
rstudio_history=filter(lambda x: x[0][0]>2013, rstudio_history)
rstudio_history=filter(lambda x: x[0][1]>7, rstudio_history)

for i in rstudio_history[:20]:
    print i[0], i[1]

#Save august 2014 Rstudio history commands
FILE_OUT1=open("RSTUDIO_FOLDER/AUGUST.2014.RSTUDIO", "w")
for i in rstudio_history:
    FILE_OUT1.write("\n"+i[1])
FILE_OUT1.close()
