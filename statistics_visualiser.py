from subprocess import Popen, PIPE
from collections import defaultdict
import matplotlib.pyplot as plt

def console(args):
    return Popen(args.split(), stdout=PIPE).communicate()[0].decode("utf-8")

def main():
    result = set()
    stdout = console("find ./cache/models -type f").split("\n")
    stdout.pop()
    number_of_models_for_id = defaultdict(int)
    number_of_models_with_rmsd = defaultdict(int)
    for model in stdout:
        model = model.split("/")
        number_of_models_for_id[model[3]] += 1
        if len(model[4]) == 3:
            continue
            #rmsd += "0"
        rmsd = float(model[4])
        number_of_models_with_rmsd[rmsd] += 1
        result.add((model[3], model[5].split(".")[0]))
    keys = []
    for el in number_of_models_with_rmsd:
        keys.append(el)
    vals = []
    for key in sorted(keys):
        #print("rmsd:", key, number_of_models_with_rmsd[key])
        vals.append(number_of_models_with_rmsd[key])

    #plt.plot(sorted(keys), vals)
    #plt.show()

    min_ = 123123123
    max_ = 0
    sum = 0
    k = 0

    number_of_models = defaultdict(int)

    for key in number_of_models_for_id:
        number_of_models[key] += number_of_models_for_id[key]
        min_ = min(min_, number_of_models_for_id[key])
        max_ = max(max_, number_of_models_for_id[key])
        sum += number_of_models_for_id[key]
        k += 1

    vals = []
    for key in number_of_models:
        vals.append(number_of_models[key])

    #plt.plot(sorted(vals))
    #plt.show()

    print("sum models:", sum)
    print("targets:", k)
    print("min models per target:", min_)
    print("max models per target:", max_)
    print("avg models per target:", sum / k)



    return result

if __name__ == "__main__":
    main()