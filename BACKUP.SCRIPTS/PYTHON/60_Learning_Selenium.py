#LEARNING selenium - REFER TO 61.py for a more in depth look at functionality

from selenium import webdriver
from selenium.webdriver.common.keys import Keys #For keys in keyboars (RETURN, F1, ALT..)

#Create the browser instance - launches it
driver=webdriver.Firefox()

#Navigate given by url to browser instance
driver.get("https://tcga-data.nci.nih.gov/tcga/dataAccessDownload.htm")
#print driver.page_source #Could print web source of some pages

#Confirm that page has title with word "Python" in it
#assert "Python" in driver.title

#Finding elements 
#elem=driver.find_element_by_name("q")
#print "#", elem.text, "#", elem.tag_name, "#", elem.location, "#", elem.parent, "#", elem.id

#Now we would send keys as if we were in the page (virtual press keys)
#elem.send_keys("selenium")
#elem.send_keys(Keys.RETURN)

#Confirm that we have reached the page we want to
#assert "Google" in driver.title

#Now you can close the browser
driver.close()

