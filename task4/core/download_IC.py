import requests
import os
import wget

from zipfile import ZipFile
url = 'http://icecube.wisc.edu/data-releases/20210126_PS-IC40-IC86_VII.zip'

wget.download(url, os.getcwd() + '/data.zip')

r = requests.get(url, allow_redirects=True)


with ZipFile(os.getcwd() + '/data.zip', 'r') as zipObj:
   # Extract all the contents of zip file in current directory
   zipObj.extractall(os.getcwd())

