"""
Submit a data query to SIMBAD for a source

http://simbad.u-strasbg.fr/simbad/sim-id?Ident=3c274&output.format=ASCII
"""

from urllib2 import urlopen
import re

def build_SIMBAD_query(source):
    url = 'http://simbad.u-strasbg.fr/simbad/'
    url += 'sim-id?Ident='
    url +=  source
    url += '&output.format=ASCII'
    return url
    
def query_SIMBAD(url):
    response = urlopen(url)
    return response

def parse_SIMBAD_response(response,
                        catalogs=['JVAS','NRAO','3C','4C','CTA','PKS','NAME']):
  html_data = response.readlines()
  in_identifiers = False
  found = {}
  for line in html_data:
    if re.search('FK5',line):
      coord_str = line.split(':')[1].strip()
      continue
    if line[:11] == "Identifiers":
      in_identifiers = True
      continue
    if in_identifiers:
      columns = line.strip().split("  ")
      clean = False
      while not clean:
        try:
          columns.remove('')
        except ValueError:
          clean=True
      for c in columns:
        column = c.strip()
        for cat in catalogs:
          if cat == column.split()[0]:
            if not found.has_key(cat):
              found[cat] = []
            found [cat].append(column[len(cat):])
    if line[:8] == "Bibcodes":
      in_identifiers = False
      break
  return found, coord_str

if __name__ == "__main__":
  source = raw_input('Enter source name:')
  url = build_SIMBAD_query(source)
  response_obj = query_SIMBAD(url)
  found, coords = parse_SIMBAD_response(response_obj)
  print found
  print coords

