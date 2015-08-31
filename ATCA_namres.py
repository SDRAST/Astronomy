"""
ATCA Source name resolver

This can be used to get the preferred source name.  The SIMBAD query will give
all the source names.

http://atcacode.googlecode.com/svn/wiki/AtnfServerSideLibrary.wiki

The mapping of variable name to string will do for now but a more general
solution is to trap the exception and then extract each variable name from that
and map it to a string, until there are no more to handle.
"""
from urllib2 import urlopen

name = 'name'
ra = 'ra'
dec='dec'
epoch = 'epoch'
position = 'position'
resolver = 'resolver'
resolvedName = 'resolvedName'

def build_URL(name, noplanets=True, noatcadb=False, nosesame=False):
  url = 'http://www.narrabri.atnf.csiro.au/cgi-bin/obstools/atnfLibrary/sourcequery.pl'
  url += '?name=' + name
  if noplanets:
    url += '&noplanets=true'
  if noatcadb:
    url += '&noatcadb=true'
  if nosesame:
    url += '&nosesame=true'
  return url

def source_query(source):
  url = build_URL(source)
  response_obj = urlopen(url)
  response = response_obj.read().strip()
  return eval(response)

if __name__ == "__main__":
  print source_query('3C274')
