#!/usr/bin/env python


modetypes = ['S','T','R']

modes = {}
for mode in modetypes:
    with open('data/modes_' + mode + '.eigen','r') as f:
        next(f)
        modes[mode] = []
        for line in f:
            line = ' '.join(line.split())
            if int(line.split(' ')[0]) == 0:
                line = line.split(' ')[4]
                modes[mode].append(1000./float(line))

for mode in modes['T']:
    print(mode)
