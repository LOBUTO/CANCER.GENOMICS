#Learning itertools.takewhile()

import itertools

X="---989SERTS9-"

Y=0
Z=[i for i in itertools.takewhile(lambda x: x=="-" or x.isdigit()==True, X)] #used to stop list comprehensions, takewhile will 
                                                                            # do what is inside the parenthesis while it is TRUE
                                                                            

print Z
print len(Z)

