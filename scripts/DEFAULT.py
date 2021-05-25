from .utilities import load, load_tar
import tarfile

class DEFAULT():
    NAME = "DEFAULT"

    CASES = ["C", "H", "random_R1", "random_R5"]

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            cases=CASES, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for case in cases:
            code = f"_{case}"
            infile = f"{input_path}{name}/{name}_{case}.pkl"
            print(f"{name} : {code}")

            loaded = load(infile)
            func(*loaded, outfile, code, exclude=[-1], timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             cases=CASES, timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for case in cases:
                code = f"_{case}"
                func1(outfile, out, { "time": t, "case": case }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             cases=CASES, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for case in cases:
            code = f"_{case}"
            infile = f"{input_path}{name}{extension}/{name}_{case}{extension}.tar.xz"
            print(f"{name} : {code}")

            tar = tarfile.open(infile)
            func(tar, timepoints, { "case": case }, outfile, code)
