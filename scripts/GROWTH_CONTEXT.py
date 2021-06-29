from .utilities import load, load_tar
import tarfile

class GROWTH_CONTEXT():
    NAME = "GROWTH_CONTEXT"

    CONTEXTS = [
        ('C', '', [-1]),
        ('CH', 'X', [-1, 4])
    ]

    POPULATIONS = ["X", "A", "B", "C", "XA", "XB", "XC", "AB", "AC", "BC",
        "XAB", "XAC", "XBC", "ABC", "XABC"]

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            contexts=CONTEXTS, populations=POPULATIONS, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for context, suffix, exclude in contexts:
            for pop in populations:
                code = f"_{context}{suffix}_{pop}"
                infile = f"{input_path}{name}/{name}_{context}_{pop}.pkl"
                print(f"{name} : {code}")

                loaded = load(infile)
                func(*loaded, outfile, code, exclude=exclude, timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             contexts=CONTEXTS, populations=POPULATIONS, timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for context, suffix, exclude in contexts:
                for pop in populations:
                    code = f"_{context}{suffix}_{pop}"
                    func1(outfile, out, { "time": t, "context": context + suffix, "pops": pop }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             contexts=CONTEXTS, populations=POPULATIONS, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for context, suffix, exclude in contexts:
            for pop in populations:
                code = f"_{context}_{pop}"
                infile = f"{input_path}{name}{extension}/{name}_{context}_{pop}{extension}.tar.xz"
                print(f"{name} : {code}")

                tar = tarfile.open(infile)
                func(tar, timepoints, { "context": context + suffix, "pops": pop }, outfile, code)
