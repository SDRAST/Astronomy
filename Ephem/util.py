def remove_item(thislist, item):
  """
  Remove item from list without complaining if absent

  @param thislist : list of items
  @param item : item which may or may not be in the list

  @return: reduced list
  """
  try:
    thislist.remove(item)
  except ValueError:
    pass
  return thislist
