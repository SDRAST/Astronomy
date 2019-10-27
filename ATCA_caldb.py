"""
Query ATCA calibrator database

http://www.narrabri.atnf.csiro.au/calibrators/calibrator_database_search.html?position=12:30:49,12:23:28&radius=1

Because it takes a while to read the database, try this:
http://stackoverflow.com/questions/28974634/urllib2-url-open-delay-between-opening-webpage-and-getting-html-data

That did not work.
"""
import urllib2
from urllib2 import urlopen
import cookielib
import re
from HTMLParser import HTMLParser
from time import sleep


class lineParser(HTMLParser):
    def __init__(self):
      HTMLParser.__init__(self)
      self.found = {}

    def handle_starttag(self, tag, attrs):
      if tag.lower() == 'h3':
        self.in_h3 = True
      else:
        pass

    def handle_endtag(self, tag):
      if tag.lower() == 'h3':
        self.in_h3 = False

    def handle_data(self, data):
      if data == "Search Results":
        self.check_result = True

def build_ATCA_caldb_query(RA='12:30:49', dec='12:23:28', radius=10):
    url = 'http://www.narrabri.atnf.csiro.au/calibrators/'
    url += 'calibrator_database_search.html?position='
    url +=  RA+','+dec
    url += '&radius='+str(radius)
    return url
    
def query_ATCA_caldb(url):
    response = urlopen(url)
    #sleep(10)
    return response

def parse_ATCA_caldb_response(response):
  html_data = response.readlines()
  parser = lineParser(catalogs=catalogs)
  for line in html_data:
    # less work for the parser to process only relevant lines
    if re.search('<A',line):
      print line
      parser.feed(line)
  return parser.found

if __name__ == "__main__":
  url = build_ATCA_caldb_query()
  cj = cookielib.CookieJar()
  opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
  req = urllib2.Request(url)
  req.add_header('User-Agent','Mozilla/5.0')
  resp = opener.open(req)
  sleep(10)
  htmltext = resp.read()
  
  response_obj = query_ATCA_caldb(url)
  found = parse_ATCA_caldb_response(response_obj)



