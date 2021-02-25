import subprocess

subprocess.run(["make","binoinv_check"])

# Starting number of samples
N = 1.0
P = [2**(-n) for n in range(1,11)]
#  m = 5
#  P = [i*2**(-m) for i in range(1,2**m)]
#  P = [1 - 2**-4]
#  P = [0.5]

results_file = "results/{prob}/binoinv_check.txt"

for p in P:
    #  if p > 0.5:
        #  continue
    output_file = results_file.format(prob=p)
    subprocess.run(["mkdir","-p", output_file.rsplit('/',1)[0]])
    subprocess.run(["./binoinv_check", str(N), str(p), output_file])
