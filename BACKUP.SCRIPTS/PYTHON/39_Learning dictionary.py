#Learning dict

D={"GA":["GAA","GAE"],"HE":["HEA","HEE","HEI"]}
print D
print D.values()

for i,j in D.iteritems():
    if "HEI" in j:
        print i

print D.has_key("FE")

l=[ [1, 'A'], [1, 'B'], [2, 'C'] ]
d={}
for key, val in l:
    d.setdefault(key, []).append(val) #d.setdefault(key,[]) returns value of key if it has it or places the next arquement
                                        #in its place (in this case an empty list) if it is not present. If there are values
                                        #such as a list, it returns the list. "val" get appended to dist list in the key values
                                        #and the dictionary is automatically updated 

print d

d[1].append('D')
print d
print d.setdefault(1,[])