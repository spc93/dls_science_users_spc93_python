import requests
#import json

#doi='10.1038/nature10414'
doi='10.1107/S2053273316000863'
url = "http://dx.doi.org/" + doi
#rapa = requests.get(url,headers={'accept':'text/x-bibliography; style=apa'})
#print rapa.text
rturt = requests.get(url,headers={'accept':'text/turtle'})
print rturt.text
#rjason = requests.get(url,headers={'accept':'application/vnd.citationstyles.csl+json'})
#print rjason.text

"""
Simple client for CrossRef DOI content negotiation
http://www.crosscite.org/cn/
After installation you can do 
> chmod +x crossRefClient.py
in order to make the file executable.

r = requests.get("http://dx.doi.org/10.1038/nature10414",headers={'accept':'text/x-bibliography; style=apa'})
print r.text

print s.doi2apa('10.1038/nature10414')


import requests
import json

#class CrossRefClient(object):

#    def __init__(self, accept='text/x-bibliography; style=apa', timeout=3):
#        """
#        # Defaults to APA biblio style
#        
#        # Usage:
#        s = CrossRefClient()
#        print s.doi2apa("10.1038/nature10414")        
#        self.headers = {'accept': accept}
#        self.timeout = timeout
#        #
#
#    def query(self, doi, q={}):
#        if doi.startswith("http://"):
#            url = doi
#        else:
#            url = "http://dx.doi.org/" + doi
#        
#        r = requests.get(url, headers = self.headers)
#        return r#
#
#
#    def doi2apa(self, doi):
#        self.headers['accept'] = 'text/x-bibliography; style=apa'
#        return self.query(doi).text
#        
#    def doi2turtle(self, doi):
#        self.headers['accept'] = 'text/turtle' 
#        return self.query(doi).text
#        
#    def doi2json(self, doi):
#        self.headers['accept'] = 'application/vnd.citationstyles.csl+json' 
#        return self.query(doi).json()
#
#'''




