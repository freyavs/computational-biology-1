

lat = []
lon = []

for i in range(17):
    f = open(f"outputs/eiwit_{i}.txt", "r")
    words = f.read()
    dd = dict()
    words = words.split('\n')
    for word in words:
        if word == "LAT": lat.append(i)
        if word == "LONG": lon.append(i)

print(lat)
print(len(lat))
print(lon)
print(len(lon))
exit()



f = open("output.txt", "r")
words = f.read()
d = dict()
words = words.split('\n')

for word in words:
    d[word] = d.get(word,0) + 1

res = []
for key in d:
    res.append( (key, d[key]) )
res = sorted(res, key= lambda tup: tup[1])

ff = open("most_common_alles.txt", "w+")
for r in res:
    ff.write(f'{r}\n')

f.close()
ff.close()

