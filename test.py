import numpy as np
from pyldpc import make_ldpc, encode, decode, get_message, coding_matrix
msg =  input()
snr = 20

binary = [int(format(ord(i), '08b')) for i in msg]
new_bin = []
for bin in binary:
    if len(str(bin)) < 8:
        bin = (8-len(str(bin))) * ["0"] + list(str(bin))
    for letter in ''.join(bin):
        new_bin.append(int(letter))

n = len(new_bin) * 2
d_v = len(new_bin) - 1
d_c = len(new_bin)

res = np.array(new_bin)
H, G = make_ldpc(n, d_v, d_c)
def readCsv():
    H = []
    f = open("matica.csv", "r")
    for line in f.readlines():
        tmp = [x for x in line.strip().split(",")]
        H.append(tmp)
    f.close()
    return H


H = np.array(readCsv()).astype(int)
G = coding_matrix(H)

y = encode(G, res, snr) # KODUJU ASI STEJNE JAKO TOTO
#print(y[int(len(y)/2):])
#print(y)

#print(''.join([str(x) for x in y]))

y = np.array([0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0 ,0 ,1 ,0, 1, 0, 1, 0, 1, 0 ,1 ,1 ,1 ,1 ,1 ,0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0])

print(y)

#if y[-1] == 1:
#    y[-1] =0 
#elif y[-1] == 0:
#    y[-1] = 1

d = decode(H, y, snr)


# overit si, ze spravne generuju - pouzit moji matici H k zakodovani a potom dekodovani

x = get_message(G, d)


def inverse(y):
    y = np.array([0 if z == 1 else 1 for z in y])
    return y

x = inverse(x)
print(x)
print(res)