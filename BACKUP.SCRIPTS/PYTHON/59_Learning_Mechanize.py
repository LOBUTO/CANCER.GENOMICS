#MECHANIZE IS NIFTY BUT UNFORTUNATELY HAVE TO GIVE UP ON IT SINCE IT CANNOT HANDLE JAVA DROP DOWN MENUS

import re
import mechanize

#Create browser object
BR=mechanize.Browser()

#Ignore robots - Important to be called before opening the page
BR.set_handle_robots(False)

BR.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
BR.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

#Open page
PAGE=BR.open("http://go.princeton.edu/cgi-bin/GOTermFinder")

#Get text of page
print PAGE.read()

print "################"
#List forms present on page
for form in BR.forms():
    print form.name #Apparently on this page there is only one form.name called "form"
    print 5
    print form #These are all the forms that are in the page that contain different controls (text input area, upload files, radiocontrol, 
                #checkbox control, submit control ...)

#To select a form when form has a name: - IMPORTANT (as stated below)
BR.select_form(nr=0) #there is only one form

#To select unnamed forms - We need to have selected a form in order to be able to select the controls in it later on

#Iterate through controls in form:
for control in BR.form.controls:
    print control
    print "type=%s, name=%s" %(control.type, control.name)


#Call specific controls - use their names
UPLOADFILE1=BR.form.find_control("uploadFile")
ONTOLOGY=BR.form.find_control("aspect")
UPLOADFILE2=BR.form.find_control("uploadBackground")
UPLOADFILE3=BR.form.find_control("uploadAnnotation")
ANNOTATION=BR.form.find_control("annoFileSelectIndex")

print 7
#Check what values can be selected in the control we called:
print UPLOADFILE1.name
print UPLOADFILE1.value #Equivalent
print BR[UPLOADFILE1.name] #Equivalent
print UPLOADFILE1.type

#Now that we know what values can be selected we can proceed to change it

#By adding files - If there are multiple fields where to add the file make sure to add the last argument specifying which control.name to add to 
BR.form.add_file(open("/Users/jzamalloa/Desktop/GO-TermFinder-0.86/examples/genes.txt"), "text/plain", 
                "/Users/jzamalloa/Desktop/GO-TermFinder-0.86/examples/genes.txt", name="uploadFile")

BR.form.add_file(open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/test"), "text/plain",
                 "/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/test", name="uploadBackground")

BR.form.add_file(open("/Users/jzamalloa/Desktop/GO-TermFinder-0.86/examples/gene_association.sgd"), "text/plain",
                 "/Users/jzamalloa/Desktop/GO-TermFinder-0.86/examples/gene_association.sgd", name="uploadAnnotation")
print UPLOADFILE1
print UPLOADFILE2
print UPLOADFILE3

#By selecting radio buttons
print ONTOLOGY.value
ONTOLOGY.value=["F"]
print ONTOLOGY

#Select from form
form.find_control("annoFileSelectIndex").readonly=False
print ANNOTATION.name
print ANNOTATION.value
print ANNOTATION.type

#When satisfied with form then submit it
RESPONSE=BR.submit(name="searchGOTermsButton")
print RESPONSE.read()



