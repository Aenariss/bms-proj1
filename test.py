import numpy as np
from pyldpc import make_ldpc, encode, decode, get_message, coding_matrix
msg =  input()
snr = 20

binary = [int(format(ord(i), '08b')) for i in msg]
new_bin = []
for bin in binary:
    for letter in str(bin):
        new_bin.append(int(letter))

n = len(new_bin) * 2
d_v = len(new_bin) - 1
d_c = len(new_bin)

res = np.array(new_bin)
#H, G = make_ldpc(n, d_v, d_c)

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

y = encode(G, res, snr)
#print(''.join([str(x) for x in y]))

#y = np.array([1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0 ,0 ,1 ,0, 1, 0, 1, 0, 1, 0 ,1 ,1 ,1 ,1 ,1 ,0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0])

############################
from BeliefPropagation import BeliefPropagation, TannerGraph, bsc_llr

#model = bsc_llr(0.3)
#tg = TannerGraph.from_biadjacency_matrix(H, channel_model=model)
#bp = BeliefPropagation(tg, H, max_iter=10)
#estimate, llr, decode_success = bp.decode(y)

#print("EST:", estimate[len(res):])


# NEUMIM OPRAVIT CHYBY
d = decode(H, y, snr)


# overit si, ze spravne generuju - pouzit moji matici H k zakodovani a potom dekodovani

x = get_message(G, d)

print("NEW:", x)