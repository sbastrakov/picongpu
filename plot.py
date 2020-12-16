import matplotlib.pyplot as plt
from math import sqrt

def read_file(name):
    f = open(name, "r")
    values = []
    lines = f.readlines()
    for line in lines:
        a, b, c = line.split()
        values.append([float(a), float(b), float(c)])
    return values

def norm(value):
    return sqrt(value[0]**2 + value[1]**2 + value[2]**2)

def get_max_norm(values):
    result = 0.0
    for value in values:
        result = max(result, norm(value))
    return result

def get_error(values, reference_values):
    max_norm = get_max_norm(reference_values)
    result = []
    for i in range(len(values)):
        first = values[i]
        second = reference_values[i]
        diff = [first[0] - second[0], first[1] - second[1], first[2] - second[2]]
        value = norm(diff) / max_norm
        result.append(value)
    return result

path_cpml_10 = "Taflove_fix3_CPML_10/simOutput/"
path_upml_10 = "Taflove_UPML_10/simOutput/"
path_reference = "Taflove_reference/simOutput/"

a_cpml_10 = read_file(path_cpml_10 + "pointA.txt")
a_upml_10 = read_file(path_upml_10 + "pointA.txt")
reference_a = read_file(path_reference + "pointA.txt")

b_cpml_10 = read_file(path_cpml_10 + "pointB.txt")
b_upml_10 = read_file(path_upml_10 + "pointB.txt")
reference_b = read_file(path_reference + "pointB.txt")

begin = 0
end = 1000

error_a_cpml_10 = get_error(a_cpml_10, reference_a)[begin:end]
error_a_upml_10 = get_error(a_upml_10, reference_a)[begin:end]
error_b_cpml_10 = get_error(b_cpml_10, reference_b)[begin:end]
error_b_upml_10 = get_error(b_upml_10, reference_b)[begin:end]

x = range(begin, end)
plt.plot(x, error_a_cpml_10, label="CPML 10 cells point A")
plt.plot(x, error_b_cpml_10, label="CPML 10 cells point B")
plt.plot(x, error_a_upml_10, label="UPML 10 cells point A")
plt.plot(x, error_b_upml_10, label="UPML 10 cells point B")
plt.xlabel('time step')
plt.ylabel('Relative error')
plt.yscale('log')
axes = plt.gca()
axes.set_ylim([1e-6,1e-2])
plt.legend()
plt.show()
