import numpy as np
from pyldpc import make_ldpc, encode, decode, get_message
msg = input()
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
print(res.shape)
print(res)
H, G = make_ldpc(n, d_v, d_c)

y = encode(G, res, snr)
d = decode(H, y, snr)
x = get_message(G, d)
print(x)