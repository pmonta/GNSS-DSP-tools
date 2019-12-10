def parse_list_ranges(s,sep='-'):
  r = []
  x = s.split(',')
  for y in x:
    z = y.split(sep)
    if len(z)==1:
      r += [int(z[0])]
    else:
      r += range(int(z[0]),int(z[1])+1)
  return list(r)

def parse_list_floats(s):
  x = s.split(',')
  return list(map(float, x))
