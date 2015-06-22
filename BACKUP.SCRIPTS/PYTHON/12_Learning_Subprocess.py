#Learning about subprocess

import subprocess

ls_output=subprocess.check_output(['ls',"-l"]) #calls command line commands, first in the LIST is the executable and the rest
                                                #are its command line arguements (-l lists permissions and size)

print ls_output #as usual needs print to print it out

ls_output


subprocess.call('ls| wc -l', shell=True) #no need for print, it will call out the shell and display the output in the console
                                        #Here is a STRING becuase we want the whole of the command to be interpreted


#subprocess.call("ls -l") #When "shell=False" won't work because there is no command called "ls -l" that can be called by
                        #the subprocess

print 6
subprocess.call(["ls", "-l"]) #Even though "shell=False" it works because it is on list format as above, calls command and then its arguements
                            #again, no need for "print"

print 7

subprocess.call("ls -l", shell=True) #again, works because shell=True, now we can use a string to call it and it will
                                    #execute as best it understands it

print 8
from subprocess import Popen

proc=Popen(["ls", "-l"])

