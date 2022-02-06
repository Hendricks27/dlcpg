import os
import sys
import matplotlib



cr = 9.575

true_re = []
pred_re = []

diff = []
for f in os.listdir("./"):

    if not f.startswith("frac"):
        continue

    if not f.endswith("test_fd42794202c5751aa48412580bf8c3cb.csv"):
        continue

    header = False
    for l in open(f):
        if not header:
            header = True
            continue
        l = l.strip().split(",")

        tre = float(l[3]) / cr
        pre = float(l[-1])

        true_re.append(tre)
        pred_re.append(pre)
        diff.append(tre - pre)


diff = list(sorted(diff, key=abs))

total_len = len(diff)

for i in [50, 75, 90, 95, 96, 97, 98, 99, 100]:
    l = int(i/100*total_len)

    newd = diff[:l]
    newd_pos = list(filter(lambda x: x >= 0, newd))
    newd_neg = list(filter(lambda x: x <  0, newd))

    mae = (sum(newd_pos) - sum(newd_neg)) / l

    p = [i] + list(map(lambda x: "%0.3f" % x, [mae, newd_neg[-1], newd_pos[-1]]))
    print("\t".join(map(str, p)))


print()
print(sum(true_re) / total_len)
print(total_len)
print(min(true_re), max(true_re))
print(min(pred_re), max(pred_re))
