def parse_list_ranges(s):
  r = []
  x = s.split(',')
  for y in x:
    z = y.split('-')
    if len(z)==1:
      r += [int(z[0])]
    else:
      r += range(int(z[0]),int(z[1])+1)
  return r

def parse_list_floats(s):
  x = s.split(',')
  return map(float, x)
