import os
import sys
import matplotlib



cr = 13.675
cr = 8.9




for f in os.listdir("./"):

    if not f.startswith("d0196_fra"):
        continue

    if "test" not in f:
        continue


    #if not f.endswith("test_fd42794202c5751aa48412580bf8c3cb.csv"):
    #    continue
    diff = []
    true_re = []
    pred_re = []

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

    print(f)
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
    print("\n---------------\n")
