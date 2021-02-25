import subprocess

# Compile with different compilers
subprocess.run(["make","timing_compilers"])
subprocess.run(["make","-C","../../../poissinv/test/Serial","timing_compilers"])

# Starting number of samples
N = 1.0
P = [2**(-n) for n in range(1,2)]

compilers = ["gcc","clang","icc"]

results_file1 = "results/{prob}/binoinv_timing_{compiler}.txt"
results_file2 = "results/{prob}/poissinv_timing_{compiler}.txt"

# Iterate over three different compilers
for compiler in compilers:
    print("\n-----------------------------------")
    print("Compiler: "+compiler.upper())
    print("-----------------------------------\n")
    for p in P:
        output_file = results_file1.format(prob=p, compiler=compiler)
        subprocess.run(["mkdir","-p", output_file.rsplit('/',1)[0]])
        subprocess.run(["./binoinv_timing_"+compiler, str(N), str(p), output_file])

        output_file = results_file2.format(prob=p, compiler=compiler)
        subprocess.run(["../../../poissinv/test/Serial/./poissinv_timing_"+compiler, str(N*p), output_file])
