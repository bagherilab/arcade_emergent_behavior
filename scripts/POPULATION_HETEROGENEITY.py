from .utilities import load, load_tar
import tarfile

class POPULATION_HETEROGENEITY():
    NAME = "POPULATION_HETEROGENEITY"

    CONTEXTS = [
        ('C', '', [-1]),
        ('CH', 'X', [-1, 4])
    ]

    POPULATIONS = ["X", "A", "B", "C", "XA", "XB", "XC", "AB", "AC", "BC",
        "XAB", "XAC", "XBC", "ABC", "XABC"]

    CANCER_HETEROGENEITY = ['00', '10', '20', '30', '40', '50']

    TISSUE_HETEROGENEITY = ['00', '10', '20', '30', '40', '50']

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            contexts=CONTEXTS, populations=POPULATIONS,
            cancer_heterogeneity=CANCER_HETEROGENEITY, tissue_heterogeneity=TISSUE_HETEROGENEITY,
            timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for context, suffix, exclude in contexts:
            for pop in populations:
                for cancer_het in cancer_heterogeneity:
                    for tissue_het in tissue_heterogeneity:
                        code = f"_{context}{suffix}_{pop}_{cancer_het}_{tissue_het}"

                        if context == "C":
                            if tissue_het != "00":
                                continue
                            infile = f"{input_path}{name}/{name}_{context}_{pop}_{cancer_het}.pkl"
                        else:
                            infile = f"{input_path}{name}/{name}_{context}_{pop}_{cancer_het}_{tissue_het}.pkl"

                        print(f"{name} : {code}")

                        loaded = load(infile)
                        func(*loaded, outfile, code, exclude=exclude, timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             contexts=CONTEXTS, populations=POPULATIONS,
             cancer_heterogeneity=CANCER_HETEROGENEITY, tissue_heterogeneity=TISSUE_HETEROGENEITY,
             timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for context, suffix, exclude in contexts:
                for pop in populations:
                    for cancer_het in cancer_heterogeneity:
                        for tissue_het in tissue_heterogeneity:
                            if context == "C" and tissue_het != "00":
                                continue

                            code = f"_{context}{suffix}_{pop}_{cancer_het}_{tissue_het}"
                            func1(outfile, out, { "time": t, "context": context + suffix, "pops": pop, "chet": cancer_het, "thet": tissue_het }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             contexts=CONTEXTS, populations=POPULATIONS,
             cancer_heterogeneity=CANCER_HETEROGENEITY, tissue_heterogeneity=TISSUE_HETEROGENEITY,
             timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for context, suffix, exclude in contexts:
            for pop in populations:
                for cancer_het in cancer_heterogeneity:
                    for tissue_het in tissue_heterogeneity:
                        code = f"_{context}_{pop}_{cancer_het}_{tissue_het}"

                        if context == "C":
                            if tissue_het != "00":
                                continue
                            infile = f"{input_path}{name}{extension}/{name}_{context}_{pop}_{cancer_het}{extension}.tar.xz"
                        else:
                            infile = f"{input_path}{name}{extension}/{name}_{context}_{pop}_{cancer_het}_{tissue_het}{extension}.tar.xz"

                        print(f"{name} : {code}")

                        tar = tarfile.open(infile)
                        func(tar, timepoints, { "context": context + suffix, "pops": pop, "chet": cancer_het, "thet": tissue_het }, outfile, code)
